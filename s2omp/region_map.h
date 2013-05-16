/*
 * region_map.h
 *
 *  Created on: Jul 24, 2012
 *      Author: cbmorrison
 */

#ifndef REGION_MAP_H_
#define REGION_MAP_H_

#include <map>
#include <stdint.h>

#include "core.h"

namespace s2omp {

typedef std::map<const uint64, int16_t> region_dict;
typedef region_dict::iterator region_iterator;
typedef std::pair<region_iterator, region_iterator> region_pair;
typedef std::map<const int16_t, double> region_area;
typedef region_area::iterator region_area_iterator;

class region_map {
  // A region_map splits a bound_interface object into roughly equal-area
  // sections.  Once that is done, one can determine the region number at a
  // given point or pixel (provided that the input pixel is at a finer
  // resolution than those used to split the bound's area).

public:
  region_map();
  virtual ~region_map();

  // By default, initialization will be done by specifying a target number
  // of regions and the code will automatically determine the appropriate
  // level to use based on the area of the bound and requested number of
  // regions.  The return value will indicate the actual number of regions
  // created.
  uint16_t init(const bound_interface& bound, uint16_t n_region);

  // Alternatively, the level for regionation can be set on initialization,
  // with the return value being the number of regions.  This method is more
  // likely to diverge from the specified number of regions if the input
  // level is not appropriate.
  uint16_t init(const bound_interface& bound, uint16_t n_region, int level);

  // Given an input point, return the region containing it (which does not
  // necessarily mean that it is contained in the bound used to generate the
  // region_map).  If the input point is outside the region map, the return
  // value is -1.
  int16_t find_region(const point& p);

  // Return the region for the input pixel.  If the pixel is outside the
  // region_map or if pix.level() < region_map.level(), the return value is -1.
  int16_t find_region(const pixel& pix);

  // Return a covering for the region indicated by the input index.
  void get_covering(int16_t region_idx, pixel_vector* pixels);

  // Given a region index, return the area associated with that region.
  double get_area(int16_t region_idx);

  // Some getter methods to describe the state of the region_map.
  inline uint16_t n_region() const {
    return n_region_;
  }
  inline uint32_t level() const {
    return level_;
  }
  inline bool is_empty() const {
    return region_map_.empty();
  }
  inline void clear() {
    region_map_.clear();
    region_area_.clear();
    level_ = -1;
    n_region_ = 0;
  }

  // Return iterators for the set of region_map objects.
  inline region_iterator begin() {
    return region_map_.begin();
  }
  inline region_iterator end() {
    return region_map_.end();
  }

private:
  int find_regionation_level(const bound_interface& bound, uint16_t n_region);
  int validate_level(double bound_area, uint16_t n_region, int level);

  region_dict region_map_;
  region_area region_area_;
  int level_;
  uint16_t n_region_;
};

} // end namespace s2omp


#endif /* REGION_MAP_H_ */
