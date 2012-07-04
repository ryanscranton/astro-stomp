// Copyright 2012  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the central class for the library: the hierarchical
// pixelization that makes all of the rest of the spatial classes work.  This
// class defines all of the core Pixel operations: converting angular position
// on the sphere into pixel index, increasing and decreasing pixel resolution,
// finding neighboring pixels and so on.

#ifndef PIXEL_H_
#define PIXEL_H_

#include <stdint.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "core.h"
#include "point.h"
#include "MersenneTwister.h"

namespace s2omp {

class point;
class pixel;

typedef std::vector<pixel> pixel_vector;
typedef pixel_vector::iterator pixel_iterator;
typedef std::pair<pixel_iterator, pixel_iterator> pixel_pair;
typedef std::vector<pixel *> pixel_ptr_vector;
typedef pixel_ptr_vector::iterator pixel_ptr_iterator;

class pixel {
  // The core class for this library.  An instance of this class represents
  // a single pixel covering a particular region of the sky, with a particular
  // weight represented by a float.  Pixels can be instantiated with an
  // AngularCoordinate and resolution level or pixel indices or just
  // instantiated with no particular location.

 public:
	explicit pixel(uint64 id);
  virtual ~pixel();

  static pixel from_point(const point& p);
  static pixel from_point(const point& p, int level);

  uint64 id() const;
  int level() const;

  bool is_leaf() const;
  bool is_face() const;
  pixel parent() const;
  pixel parent(int level) const;
  void children(pixel_vector& child_pixels) const;
  void children(int level, pixel_vector& child_pixels) const;
  pixel child_begin() const;
  pixel child_begin(int level) const;
  pixel child_end() const;
  pixel child_end(int level) const;

  pixel next() const;
  pixel next_wrap() const;
  pixel prev() const;
  pixel prev_wrap() const;

  static double average_area(int level);
  double average_area() const;
  double exact_area() const;

  bool contains(point& p) const;
  bool contains(pixel& pix) const;

  point center_point() const;
  point vertex(int k) const;
  point edge(int k) const;

  pixel_vector* neighbors() const;
  pixel_vector* neighbors(int level) const;

 private:
  pixel();

  S2::S2CellId id_;
  S2::S2Cell* cell_;
};

} // end namespace s2omp

#endif /* PIXEL_H_ */
