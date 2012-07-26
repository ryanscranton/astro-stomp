/*
 * region_map.cc
 *
 *  Created on: Jul 24, 2012
 *      Author: cbmorrison
 */

#include "region_map.h"

namespace s2omp {

region_map::region_map() {
  clear();
}

region_map::~region_map() {
  clear();
}

uint16_t region_map::init(const bound_interface& bound, uint16_t n_region) {
  return init(bound, n_region, find_regionation_level(bound, n_region));
}

uint16_t region_map::init(const bound_interface& bound, uint16_t n_region,
    int level) {
  clear();

  // First, verify that our input level is valid.
  double bound_area = bound.area();
  level = validate_level(bound_area, n_region, level);

  // Find the pixels that cover this bound at the level requested.
  pixel_vector pixels;
  coverer.simple_covering(bound, level, &pixels);

  // Since coverings aren't guaranteed to give us an ordered list of pixels we
  // need to sort the covering pixels.
  sort(pixels.begin(), pixels.end());

  uint16_t created_regions = 0;
  double region_area = 0.0;
  double pixel_average_area = pixel.average_area(level);
  double target_region_area = bound_area / (1.0 * n_region);
  for (pixel_iterator iter = pixels.begin(); iter != pixels.end(); ++iter) {
    // For each covering pixel we need to know the how much area overlaps with
    // this bound.
    double contained_area = bound.contained_area(*iter);

    // This is the test for sending the pixel to the current or next region.
    // The test is fuzzy in the sense that as long as we are within the average
    // area of one pixel of the target area we complete the region. If we are
    // on the last region, just keep adding until we reach the end.
    if (region_area + contained_area < target_region_area - pixel_average_area
        || created_regions == n_region - 1) {
      region_map_.insert(std::pair<const uint64, uint16_t>(iter->id(),
          created_regions));
      region_area += contained_area;
    } else {
      // Store current region area and add this pixel to the next region,
      // reseting the region area.
      region_area_.insert(std::pair<const uint16_t, double>(created_regions,
          region_area));

      created_regions++;
      region_map_.insert(std::pair<const uint64, uint16_t>(iter->id(),
          created_regions));
      region_area = contained_area;
    }
  }
  region_area_.insert(std::pair<const uint16_t, double>(created_regions,
      region_area));
  created_regions++;

  level_ = level;
  n_region_ = created_regions;

  pixels.clear();

  return n_region_;
}

int16_t region_map::find_region(const point& p) {
  region_iterator iter = region_map_.find(p.id(level_));
  if (iter != end())
    return iter->second;
  return -1;
}

int16_t region_map::find_region(const pixel& pix) {
  if (pix.level() >= level_) {
    region_iterator iter = region_map_.find(pix.parent(level_).id());
    if (iter != end()) {
      return iter->second;
    }
  }

  return -1;
}

void region_map::region_covering(int16_t region, pixel_vector* pixels) {
  if (!pixels->empty())
    pixels->clear();

  for (region_iterator iter = begin(); iter != end(); ++iter) {
    if (iter->second == region) {
      pixels->push_back(pixel(iter->first));
    } else if (iter - second > region) {
      break;
    }
  }
}

double region_map::region_area(int16_t region) {
  region_iterator iter = region_area_.find(region);
  if (iter != end())
    return iter->second;
  return -1.0;
}

int region_map::find_regionation_level(const bound_interface& bound,
    uint16_t n_region) {
  double target_area = bound.area() / (50 * n_region);
  int level = 0;
  while ((pixel::average_area(level) > target_area) && (level < MAX_LEVEL)) {
    level++;
  }

  return level;
}

int region_map::validate_level(double bound_area, uint16_t n_region,
    int level) {
  if (level > MAX_LEVEL) {
    std::cout << "region_map::init - Requested level > MAX_LEVEL."
        << "\tExiting regionation.";
    exit(2);
  }

  if (level < 0) {
    std::cout << "region_map::init - Requested level < 0"
        << "\tExiting regionation.";
    exit(2);
  }

  // In order for regionation to generate the number of requested regions,
  // we need at least as many covering pixels as requested regions.  If that's
  // not the case, increase the level until we satisfy that requirement.
  double target_region_area = bound_area / (1.0 * n_region);
  while (target_region_area < pixel.average_area(level)) {
    level++;
    std::cout << "region_map::init - Average pixel size at requested "
        << "level is greater than expected region area.\n"
        << "\nincreasing level to " << level << ".\n";
    if (level > MAX_LEVEL) {
      std::cout << "region_map::init - Regionation level too high"
          << "\tfailing regionation.";
      exit(2);
    }
  }

  return level;
}

} // end namespace s2omp
