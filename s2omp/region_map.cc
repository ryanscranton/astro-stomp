/*
 * region_map.cc
 *
 *  Created on: Jul 24, 2012
 *      Author: cbmorrison
 */

#include "bound_interface.h"
#include "coverer.h"
#include "pixel.h"
#include "point.h"
#include "region_map.h"

namespace s2omp {

region_map::region_map() {
  clear();
}

region_map::~region_map() {
  clear();
}

uint region_map::init(const bound_interface& bound, uint n_region) {
  return init(bound, n_region, find_regionation_level(bound, n_region));
}

uint region_map::init(const bound_interface& bound, uint n_region, int level) {
  clear();

  // First, verify that our input level is valid.
  double bound_area = bound.area();
  level = validate_level(bound_area, n_region, level);

  // Find the pixels that cover this bound at the level requested.  The covering
  // pixels will be returned sorted.
  pixel_vector pixels;
  bound.get_simple_covering(level, &pixels);

  uint created_regions = 0;
  double region_area = 0.0;
  double pixel_average_area = pixel::average_area(level);
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
      region_map_[iter->id()] = created_regions;
      region_area += contained_area;
    } else {
      // Store current region area and add this pixel to the next region,
      // reseting the region area.
      region_area_[created_regions] = region_area;

      created_regions++;
      region_map_[iter->id()] = created_regions;
      region_area = contained_area;
    }
  }
  region_area_[created_regions] = region_area;
  created_regions++;

  level_ = level;
  n_region_ = created_regions;

  pixels.clear();

  return n_region_;
}

int region_map::find_region(const point& p) const {
  region_iterator iter = region_map_.find(p.id(level_));
  // If the point is within our region_map, return the region.
  return iter != region_map_.end() ? iter->second : INVALID_REGION_VALUE;
}

int region_map::find_region(const pixel& pix) const {
  // If the pixel is larger than our region_map level.
  if (pix.level() < level_) {
    return INVALID_REGION_VALUE;
  }

  // If the pixel is within our region_map, return the region.
  region_iterator iter = region_map_.find(pix.parent(level_).id());
  return iter != region_map_.end() ? iter->second : INVALID_REGION_VALUE;
}

int region_map::find_region(const region_map& regions, const point& p) {
  int region = !regions.is_empty() ? regions.find_region(p)
      : INVALID_REGION_VALUE;
  if (!regions.is_empty() && region == INVALID_REGION_VALUE) {
    std::cout << "s2omp::region_map::find_region - "
        << "Failed to find region for input point; bailing\n";
    exit(2);
  }

  return region;
}

int region_map::find_region(const region_map& regions, const pixel& pix) {
  int region = !regions.is_empty() ? regions.find_region(pix)
      : INVALID_REGION_VALUE;
  if (!regions.is_empty() && region == INVALID_REGION_VALUE) {
    std::cout << "s2omp::region_map::find_region - "
        << "Failed to find region for input pixel; bailing\n";
    exit(2);
  }

  return region;
}

void region_map::get_covering(int region_idx, pixel_vector* pixels) const {
  if (!pixels->empty()) {
    pixels->clear();
  }

  // If the input region is outside our range, stop now.
  if (region_idx < 0 || region_idx >= n_region_) {
    return;
  }

  for (region_iterator iter = region_map_.begin(); iter != region_map_.end(); ++iter) {
    if (iter->second == region_idx) {
      pixels->push_back(pixel(iter->first));
    }
  }
}

double region_map::get_area(int region) const {
  // If the region index is valid, return the area.
  region_area_iterator iter = region_area_.find(region);
  return iter != region_area_.end() ? iter->second : -1.0;
}

int region_map::find_regionation_level(const bound_interface& bound,
    uint n_region) const {

  // To cleanly sub-divide our bound area (i.e., each region has approximately
  // the same area), we want on order 50 pixels per region.  Increase the level
  // until we get to that resolution.
  double target_area = bound.area() / (50 * n_region);
  int level = 0;
  while ((pixel::average_area(level) > target_area) && (level < MAX_LEVEL)) {
    level++;
  }

  return level;
}

int region_map::validate_level(double bound_area, uint n_region, int level) const {
  // TODO(scranton): Refactor this so that we fail more gracefully.
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
  while (target_region_area < pixel::average_area(level)) {
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
