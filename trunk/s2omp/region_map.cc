/*
 * region_map.cc
 *
 *  Created on: Jul 24, 2012
 *      Author: cbmorrison
 */


#include "region_map.h"

namespace s2omp {

region_map::region_map() {
  level_ = -1;
  n_region_ = 0;
  region_map_.clear();
  region_area_.clear();
}

uint16_t region_map::initialize(uint16_t n_region,
    const bound_interface& bound) {
  find_regionation_level(bound, n_region);
  return initialize(n_region, level, bound);
}

uint16_t region_map::initialize(uint16_t n_region, int level,
    const bound_interface& bound) {
  region_map_.clear();
  region_area_.clear();

  // before doing anything test that the requested level is not greater than
  // MAX_LEVEL
  if (level > MAX_LEVEL) {
    std::cout << "region_map::initialize - Requested level > MAX_LEVEL."
        << "\tExiting regionation.";
    exit(2);
  }

  // to test if the regionation request is sane and keep track of how much area
  // we have regionated, we compute the total area and the target area for each
  // of the regions.
  double total_area = bound.area();
  double target_region_area = total_area/(1.0*n_region);

  // want to test that the user input is sane.
  while (target_region_area < pixel.average_area(level)) {
    level++;
    std::cout << "region_map::initialize - Average pixel size at requested "
        << "level is greater than expected region area.\n"
        << "\nincreasing level to " << level << ".\n";
    if (level > MAX_LEVEL_) {
      std::cout << "region_map::initialize - Regionation level too high"
          << "\tfailing regionation.";
      exit(2);
    }
  }
  // now that we have the regionation level we can compute the average pixel
  // area for use in testing when to switch regions.
  double average_area = pixel.average_area(level);

  // find the pixels that cover this bound at the level requested.
  pixel_vector* pixels;
  coverer.simple_covering(level, bound, pixels);

  // since coverings aren't guaranteed to give us an ordered list of pixels we
  // need to sort covering's output.
  sort(pixels->begin(), pixels->end());

  // book keeping variables
  uint16_t created_regions = 0;
  double region_area = 0.0;

  for (pixel_iterator iter = pixels->begin(); iter != pixels->end(); ++iter) {
    // for each covering pixel we need to know the how much area overlaps with
    // this bound
    double area = bound.contained_area(*iter);
    // this is the test for sending the pixel to the current or next region.
    // the test is fuzzy in the sense that as long as we are within the average
    // area of one pixel of the target area we complete the region. If we are on
    // the last region, just keep adding until we reach the end.
    if (region_area + area < target_region_area - average_area ||
        created_regions == n_region - 1) {
      region_map_.insert(std::pair<const uint64, uint16_t>(
          iter->id(), created_regions));
      region_area += area;
    } else {
      // store current region area and add this pixel to the next region,
      // reseting the region area.
      region_area_.insert(
          std::pair<const uint16_t, double>(created_regions, region_area));
      created_regions++;
      region_area = 0.0;

      region_map_.insert(std::pair<const uint64, uint16_t>(
          iter->id(), created_regions));
      region_area += area;
    }
  }
  region_area_.insert(
      std::pair<const uint16_t, double>(created_regions, region_area));
  created_regions++;

  level_ = level;
  n_region_ = created_regions;

  pixels.clear();
  delete pixels;

  return created_regions;
}

bool region_map::initialize(const region_map regions,
    bound_interface* bound) {
  if (!regions.is_initialized()) return false;

  // check the covering on the bound. If it does not contain the same number of
  // pixels as this bound then we know these bounds do not have the same area
  // (Necessary but not sufficient)
  pixel_vector* pixels;
  coverer.simple_covering(level_, bound, pixels);
  if (pixels->size() != region_map_.size()) return false;

  bool copied_regions = true;
  // Now we need to test if the pixels in the covering are the same geometry
  for (pixel_iterator iter = pixels->begin(); iter != pixels->end(); ++iter) {
    region_iterator region_iter = region_map_.find(iter->id());
    if (region_iter == end()) {
      copied_regions = false;
      break;
    }
  }

  pixels.clear();
  delete pixels;

  return copied_regions;
}

int16_t region_map::find_region(const point& p) {
  region_iterator iter = region_map_.find(p.id(level_));
  if (iter != end())
    return iter->second;
  return -1;
}

int16_t region_map::find_region(const pixel& pix) {
  int16_t region = -1;
  if (pix.level() < level_) {
    region_iterator iter = region_map_.find(pix.parent().id());
    if (iter != end())
      region = iter->second;
  } else if (pix.level() == level_) {
    region_iterator iter = region_map_.find(pix.id());
    if (iter != end())
      region = iter->second;
  } else {
    pixel_vector* pixels;
    pix.children(level_, pixels);

    for (pixel_iterator iter = pixels->begin(); iter != pixels->end(); ++iter) {
      region_iterator region_iter = region_map_.find(iter->id());
      if (region_iter == end()) {
        region = -1;
        break;
      } else if (region_iter->second != region || region != -1) {
        region = -1;
        break;
      } else if (region == -1) {
        region = region_iter->second;
      }
    }
    pixels->clear();
    delete pixels;
  }
  return region;
}

void region_map::region_covering(int16_t region, pixel_vector* pixels) {
  if (!pixels->empty()) pixels->clear();

  for (region_iterator iter = being(); iter != end(); ++iter) {
    if (iter->second == region) {
      pixels->push_back(pixel(iter->first));
    } else if (iter-second > region) {
      break;
    }
  }
}

double region_map::region_area(int16_t region) {
  region_iterator iter = region_area_.find(region);
  if (iter != end()) return iter->second;
  return -1.0;
}

int region_map::find_regionation_level(const bound_interface& bound,
    uint16_t n_region) {
  double target_region_area = bound.area()/(1.0*n_region);
  int level = 0;
  while (target_region_area < pixel.average_area(level)) {
    level++;
    if (level > MAX_LEVEL) {
      std::cout << "region_map::initialize - Regionation level too high"
          << "\tfailing regionation.";
      exit(2);
    }
  }

  return level;
}

} // end namespace s2omp
