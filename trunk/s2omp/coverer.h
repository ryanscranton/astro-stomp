/*
 * coverer.h
 *
 *  Created on: Jul 24, 2012
 *      Author: cbmorrison
 */

#ifndef COVERER_H_
#define COVERER_H_

#include <queue>
#include <set>
#include <utility>

#include "core.h"
#include "pixel.h"

namespace s2omp {

class pixel;

struct pixel_candidate {
  pixel pix;
  bool is_terminal;
  uint8_t n_children;
};

typedef std::pair<int, pixel_candidate> pixel_entry;
typedef std::vector<pixel_entry> pixel_entry_vector;
typedef pixel_entry_vector::iterator pixel_entry_iterator;

struct compare_pixels : public std::less<pixel_entry> {
  bool operator()(pixel_entry const& x, pixel_entry const& y) {
    return x.first < y.first;
  }
};

typedef std::priority_queue<pixel_entry, pixel_entry_vector,
                            compare_pixels> pixel_queue;

class coverer {

// A Coverer is a class that approximates and arbitrary bound bound object
// by pixels that "cover" the bound. The level of the pixels returned are
// between min_level - max_level as specified at construction. This class is
// what is used to create the initial set of pixels to create a pixel_union.

public:
  //Default constructor. Sets min_level to 0 and max_level to MAX_LEVEL
  coverer();
  // Constructor for explicitly stating a min and max level. For interior
  // interior coverings it is recommend that a max level be specified.
  coverer(int min_level, int max_level);
  virtual ~coverer();

  // A covering returns a set of pixels that cover the geometry specified in
  // bound, approximating the area. The returned covering area will be greater
  // than the bound area. The two versions below operate by specifying
  // a maximum number of pixels to use in the covering (default 8). A covering
  // should always return a vector of pixels that satisfy may intersect.
  // To give the user an idea of how many pixels to specify, here is the
  // documentation from S2. Note that cells are equivalent to pixels.

  // "Accuracy is measured by dividing the area of the covering by the area of
  // the original region.  The following table shows the median and worst case
  // values for this area ratio on a test case consisting of 100,000 spherical
  // caps of random size (generated using s2regioncoverer_unittest):
  //
  //   max_cells:        3      4     5     6     8    12    20   100   1000
  //   median ratio:  5.33   3.32  2.73  2.34  1.98  1.66  1.42  1.11   1.01
  //   worst case:  215518  14.41  9.72  5.26  3.91  2.75  1.92  1.20   1.02"
  bool get_covering(const bound_interface& bound, pixel_vector* pixels);
  bool get_covering(const uint32_t max_pixels, const bound_interface& bound,
      pixel_vector* pixels);

  // This covering works by specifying a target fractional area tolerance. The
  // method will continually refine the covering until either the tolerance is
  // achieved or no pixels at a level above max_level_ (that are not contained)
  // remain. The return is a boolean specifying if the tolerance was reached.
  bool get_covering(double fractional_area_tolerace,
    const bound_interface& bound, pixel_vector* pixels);

  // Interior coverings are different from the coverings above in that each
  // pixel in the covering is required to be fully contained within the bound.
  // Here, no max number of pixels is specified. The method will terminate once
  // all pixels at max_level_ that may intersect the bound are tested for
  // contains. For interior coverings, it is possible that the return pixels
  // will be empty if the max level specified has no pixels that are fully
  // contained within the bound.
  bool get_interior_covering(const bound_interface& bound,
                             pixel_vector* pixels);

  // Similar to above only this time in addition to the max level, the method
  // will also terminate upon reaching a specified tolerance. The return boolean
  // specifies if this tolerance has been met.
  bool get_interior_covering(
      double fractional_area_tolerace, const bound_interface& bound,
      pixel_vector* pixels);

  // A simple covering is a set of contiguous at a single level that cover
  // a bound. Note that if bound is a discontinuous region on the sky or has
  // large holes in the bound, it is possible that for certain levels that
  // simple covering will not return a covering for the whole area or will
  // return an empty vector of pixels.
  static void get_simple_covering(const bound_interface& bound, int level,
                                  pixel_vector* pixels);

  inline int get_min_level() {
    return min_level_;
  }
  inline int get_max_level() {
    return max_level_;
  }

  inline bool set_min_level(int level) {
    return set_min_max_level(level, max_level_);
  }
  inline bool set_max_level(int level) {
    return set_min_max_level(min_level_, level);
  }
  bool set_min_max_level(int min, int max);

private:
  bool generate_covering(const bound_interface& bound, uint32_t max_pixels,
                         bool interior, double fraction, pixel_vector* pixels);
  void get_initial_covering(const bound_interface& bound, pixel_vector* pixels);
  int score_pixel(const bound_interface& bound, pixel_candidate* pix);
  void flush_queue(pixel_vector* pixels);
  void new_candidate(const bound_interface& bound, const pixel& pix);

  int min_level_, max_level_;
  pixel_queue pix_q_;
};

} // end namespace s2omp


#endif /* COVERER_H_ */
