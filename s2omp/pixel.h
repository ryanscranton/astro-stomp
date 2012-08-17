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

#include "circle_bound.h"
#include "point.h"
#include "MersenneTwister.h"

// S2 Includes for pixel.h
#include "s2.h"
#include "s2cellid.h"
#include "s2cell.h"

namespace s2omp {

typedef std::vector<pixel> pixel_vector;
typedef pixel_vector::iterator pixel_iterator;
typedef std::pair<pixel_iterator, pixel_iterator> pixel_pair;
typedef std::vector<pixel *> pixel_ptr_vector;
typedef pixel_ptr_vector::iterator pixel_ptr_iterator;

class circle_bound;
class pixel;
class point;

typedef std::vector<pixel> pixel_vector;
typedef pixel_vector::iterator pixel_iterator;
typedef std::pair<pixel_iterator, pixel_iterator> pixel_pair;
typedef std::vector<pixel *> pixel_ptr_vector;
typedef pixel_ptr_vector::iterator pixel_ptr_iterator;

class pixel {
  // The core class for this library.  An instance of this class represents
  // a single pixel covering a particular region of the sky.  The bulk of the
  // class functionality is based on an S2::S2CellId, but it also combines some
  // the S2::S2Region-like behavior of the S2::S2Cell class.

public:
  // The default constructor requires a 64-bit integer (based on the S2
  // pixelization).
  explicit pixel(uint64 id) : id_(S2CellId(id)) {}
  virtual ~pixel();

  // Static constructors are also available to instantiate
  // from point objects (either at a specific input level or at the maximum
  // resolution).
  static pixel* from_point(const point& p);
  static pixel* from_point(const point& p, int level);

  // Accessor methods for our basic parameters.
  inline uint64 id() const {
    return id_.id();
  }
  inline int level() const {
    return id_.level();
  }

  // Simple methods for checking whether our pixel is at the top or the
  // bottom of our resolution limits (roughly 1/6th of the sky down to about
  // 0.00015 arcseconds).
  inline bool is_leaf() const {
    return id_.is_leaf();
  }
  inline bool is_face() const {
    return id_.is_face();
  }

  // Our basic methods for transversing the pixelization hierarchy.  In both
  // cases, we can either access the parent or child cells directly above or
  // below our current level or at an input level.
  pixel parent() const;
  pixel parent(int level) const;
  void children(pixel_vector* child_pixels) const;
  void children(int level, pixel_vector* child_pixels) const;

  // Alternatively, if want to avoid allocating the memory for our child
  // pixels in one go, we can iterate through them.  The syntax for this is
  // similar to the corresponding methods for S2::S2CellId:
  //
  // for (pixel c = pix.child_begin(); c != pix.child_end(); c = c.next()) {
  // ...
  //
  // As with S2::S2CellId, the iterator method should be called rather than
  // ++c.
  pixel child_begin() const;
  pixel child_begin(int level) const;
  pixel child_end() const;
  pixel child_end(int level) const;

  // Iterator methods for moving along the Hilbert curve at our pixel's current
  // level.  The *_wrap() methods allow the iterators to move across the pixel
  // index discontinuity.
  pixel next() const;
  pixel prev() const;
  pixel next_wrap() const;
  pixel prev_wrap() const;

  // Our first methods using the functionality from S2::S2Cell.  Since S2 cells
  // are only roughly equal-area, there is some small difference between the
  // average area of a cell at a given level at the exact area of a given cell.
  double average_area(int level);
  double average_area();
  double exact_area();

  // For some methods we need would like a way to calculate the level at which
  // the average area of the level is less than the input area. This function
  // will return levels less than 0 and greater than 30 so it must be tested
  // after calling if a valid level is needed.
  inline static int get_level_from_area(double area) {
    return int(ceil(log(21600.0/(PI*area))/log(4.0)));
  }

  // Return true if the input point or pixel is inside our current pixel.
  bool contains(const point& p) const;
  bool contains(const pixel& pix) const;

  // In addition, we can also check to see if we are contained by another
  // pixel.
  bool intersects(const pixel& pix) const;

  // The real work for these containment methods is actually done by
  // range_min() and range_max(), which return the pixels that bound our pixel
  // at all levels.  If a second pixel has an id between range_min() and
  // range_max(), then we know that it is contained by our pixel.
  inline pixel range_min() const;
  inline pixel range_max() const;

  // Methods for extracting the points that define the center, vertices and
  // edges of our pixel.  The pixel edges are defined by great circles, so the
  // returned points are the vectors orthogonal to a given great circle,
  // pointed interior to the pixel.  The ordering of the vertices and edges
  // are such that edge(0) runs from from vertex(0) to vertex(1).
  point center_point() const;
  point vertex(int k) const;
  point edge(int k) const;

  // For some uses of the code we would like to know both the distance to
  // the closest and farthest edge given an external point.
  double nearest_edge_distance(const point& p);
  double farthest_edge_distance(const point& p);

  // Alternatively, we can get both edge distances in a single call.  The
  // returned boolean tells us whether the distances are to pixel edges (true)
  // or pixel corners (false).
  bool edge_distances(const point& p, double& near_edge_distance,
      double& far_edge_distance);

  // Return a vector of pixels containing the neighbors of this pixel.  In the
  // second case, we return the neighboring pixels at an input level.
  void neighbors(pixel_vector* pixels) const;
  void neighbors(int level, pixel_vector* pixels) const;

  // Return a Poisson-random point (or multiple such points) from the area
  // covered by this pixel
  point get_random_point() const;
  void get_random_points(long n_points, pixel_vector* points) const;

protected:
  void set_id(uint64 id);

private:
  pixel();
  S2Cell get_cell() const;
  static point s2point_to_point(const S2Point& p);
  static S2Point point_to_s2point(const point& p);

  S2CellId id_;
};

inline bool operator==(pixel const& a, pixel const& b) {
  return a.id() == b.id();
}

inline bool operator!=(pixel const& a, pixel const& b) {
  return a.id() != b.id();
}

inline bool operator<(pixel const& a, pixel const& b) {
  return a.id() < b.id();
}

inline bool operator>(pixel const& a, pixel const& b) {
  return a.id() > b.id();
}

inline bool operator<=(pixel const& a, pixel const& b) {
  return a.id() <= b.id();
}

inline bool operator>=(pixel const& a, pixel const& b) {
  return a.id() >= b.id();
}

} // end namespace s2omp

#endif /* PIXEL_H_ */
