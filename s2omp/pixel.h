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
  pixel(uint64 id);

  static pixel* from_point(point& p);
  static pixel* from_point(point& p, int level);

  uint64 id();
  int level();

  bool is_leaf();
  bool is_face();
  pixel parent();
  pixel parent(int level);
  pixel_vector* children();
  pixel_vector* children(int level);
  pixel child_begin();
  pixel child_begin(int level);
  pixel child_end();
  pixel child_end(int level);

  pixel next();
  pixel next_wrap();
  pixel prev();
  pixel prev_wrap();

  static double average_area(int level);
  double average_area();
  double exact_area();

  bool contains(point& p);
  bool contains(pixel& pix);

  point center_point();
  point vertex(int k);
  point edge(int k);

  pixel_vector* neighbors();
  pixel_vector* neighbors(int level);

 private:
  pixel();
  virtual ~pixel();

  S2::S2CellId id_;
  S2::S2Cell* cell_;
};

} // end namespace s2omp

#endif /* PIXEL_H_ */
