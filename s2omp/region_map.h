/*
 * region_map.h
 *
 *  Created on: Jul 24, 2012
 *      Author: cbmorrison
 */

#ifndef REGION_MAP_H_
#define REGION_MAP_H_

typedef std::map<const uint64, int16_t> region_dict;
typedef region_dict::iterator region_iterator;
typedef std::pair<region_iterator, region_iterator> region_pair;
typedef std::map<const int16_t, double> region_area;

namespace s2omp {
class region_map {
  // A region_map splits a bound_interface object into roughly equal-area
  // sections.  Once that is done, one can determine the region number at a
  // given point or pixel (provided that the input pixel is at a finer
  // resolution than those used to split the bound's area).

public:
  region_map();
  virtual ~region_map();

  // There are two options for
  uint16_t init(uint16_t n_region, bound_interface* bound);
  uint16_t init(uint16_t n_region, int level, bound_interface* bound);

  int16_t find_region(const point& p);

  int16_t find_region(const pixel& pix);

  void clear();

  void region_covering(int16_t region, pixel_vector* pixels);

  // Given a region index, return the area associated with that region.
  double region_area(int16_t region);

  // Some getter methods to describe the state of the RegionMap.
  inline uint16_t n_region() const {
    return n_region_;
  }
  inline uint32_t level() const {
    return level_;
  }
  inline bool is_initialized() const {
    return !region_map_.empty();
  }

  // Return iterators for the set of RegionMap objects.
  inline region_iterator begin() {
    return region_map_.begin();
  }
  inline region_iterator end() {
    return region_map_.end();
  }

private:
  int find_regionation_level(const bound_interface& bound, uint16_t n_region);

  void regionate(const bound_interface& bound, uint16_t n_region, int level);

  bool verify_regionation(uint16_t n_region);

  region_dict region_map_;
  region_area region_area_;
  int level_;
  uint16_t n_region_;
};

} // end namespace s2omp


#endif /* REGION_MAP_H_ */
