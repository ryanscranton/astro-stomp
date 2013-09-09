// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the TreeMap class.  The core work of pair finding
// and K nearest neighbor searches is done in the TreePixel class.  However,
// due to the pixelization scheme used in STOMP, there is a maximum pixel size
// that does not span the entire sphere.  Hence, a vector of TreePixels is
// necessary to describe an arbitrary collection of points.  TreeMap manages
// that vector of TreePixels, adding them as necessary based on the input
// points.

#ifndef TREE_UNION_H_
#define TREE_UNION_H_

#include <stdint.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <algorithm>
#include "core.h"
#include "bound_interface.h"
#include "circle_bound.h"
#include "point.h"
#include "region_map.h"
#include "tree_pixel.h"

namespace s2omp {
typedef std::map<const uint64, tree_pixel *> tree_map;
typedef tree_map::const_iterator tree_map_iterator;
typedef std::pair<tree_map_iterator, tree_map_iterator> tree_map_pair;
typedef std::pair<tree_map_iterator, bool> tree_map_insert_iterator;

class tree_union: public bound_interface {
public:
  friend class nearest_neighbor_pixel;

  tree_union();
  tree_union(int level);
  tree_union(int level, int max_points);

  bool add_point(const point& p);

  // Stand-alone methods where we're only interested in the pairs and/or
  // weight within a single bound.
  long find_pairs(const annulus_bound& bound) const;
  double find_weighted_pairs(const annulus_bound& bound) const;

  // If we want to find pairs within an angular bin for a pair-based
  // correlation function, then we have a different interface that stores the
  // results in an angular_bin object that also defines the annulus around each
  // of the points in the input vector.  In this case, we store both the
  // total number of pairs (in bin->pair_counts()) and product sum of the
  // weights (in bin->pair_weight()).  The return value indicates success or
  // failure in finding all of the expected pairs.
  bool find_pairs(const point_vector& points, angular_bin* bin) const;
  bool find_pairs_with_regions(const point_vector& points,
      const region_map& regions, angular_bin* bin) const;

  // In addition to pair finding, we can also use the tree structure we've
  // built to do efficient nearest neighbor searches.  In the general case,
  // we'll be finding the k nearest neighbors of an input point.  The return
  // value is the number of nodes touched during the assemblage.
  //
  // NOTE: There is no duplication checking.  Hence, if the input point is a
  // copy of a point in the tree, then that point will be included in the
  // returned vector of points.
  //
  // The central engine of our neighbor finding uses recursion through the
  // nodes to find nearest matches.  However, we generally don't want to use
  // the method directly, so we start the method with an underscore to indicate
  // its semi-private nature.
  void _neighbor_recursion(const point& p, tree_neighbor* neighbor) const;

  // Our first preferred method, return a vector of the nearest k neighbors
  // to the input point.
  long find_k_nearest_neighbors(const point& p, int n_neighbors,
      point_vector* neighbors) const;

  // The special case where we're only interested in the nearest matching point.
  long find_nearest_neighbor(const point& p, point* neighbor) const;

  // In some cases, we're only interested in the distance to the kth nearest
  // neighbor.  The return value will be the angular distance in degrees.
  double k_nearest_neighbor_distance(const point& p, int n_neighbors,
      long& nodes_visited) const;

  // Or in the distance to the nearest neighbor.
  double nearest_neighbor_distance(const point& p, long& nodes_visited) const;

  // Alternatively, we could be less interested in the nearest neighbor and
  // more interested in finding a direct match to our input point.  The
  // difference is subtle, but whereas NearestNeighbor will always return a
  // point from our tree, ClosestMatch has an angular threshold, beyond which
  // we're not interested in the nearest neighbor because it's not a match to
  // our input point.  The returned boolean indicates whether the returned
  // point is within the specified radius and the maximum distance is
  // in degrees.  As with nearest neighbors, the preferred interface is
  // "closest_match" and not "_match_recursion".
  void _match_recursion(const point& p, tree_neighbor* neighbor) const;
  bool
      closest_match(const point& p, double max_angular_distance, point& match) const;

  long n_points() const {
    return point_count_;
  }
  double weight() const {
    return weight_;
  }

  // A variation on the above method, returns the number of points associated
  // with the tree_union that are also contained in the input pixel.
  long n_points(const pixel& pix) const;
  double weight(const pixel& pix) const;

  inline int level() const {
    return level_;
  }
  inline int pixel_capacity() const {
    return maximum_points_;
  }

  // If we want to extract a copy of all of the points that have been added
  // to this pixel, this method allows for that.
  void points(point_vector* points) const;

  // And an associated method that will extract a copy of the points associated
  // with an input pixel.
  void points(const pixel& pix, point_vector* points) const;

  // API from pixelized_bound_interface.h
  inline virtual bool is_empty() const {
    return !tree_map_.empty();
  }
  inline virtual long size() const {
    return tree_map_.size();
  }
  virtual void clear();
  virtual double area() const;

  virtual bool contains(const point& p) const;
  virtual bool contains(const pixel& pix) const;

  virtual double contained_area(const pixel& pix) const;
  virtual bool may_intersect(const pixel& pix) const;

  virtual circle_bound get_bound() const;
  virtual point get_center() const;

  virtual void get_covering(pixel_vector* pixels) const;
  virtual void get_size_covering(const long max_pixels, pixel_vector* pixels) const;
  virtual void get_area_covering(double fractional_area_tolerance,
      pixel_vector* pixels) const;
  virtual void get_interior_covering(int max_level, pixel_vector* pixels) const;
  virtual void get_simple_covering(int level, pixel_vector* pixels) const;
  virtual void get_center_covering(int level, pixel_vector* pixels) const;

private:
  tree_map_iterator add_node(uint64 id);
  void get_node_level_pixels(const pixel& pix, pixel_vector* pixels) const;
  void initialize_bound();
  void calculate_area();

  tree_map tree_map_;
  pixel_set nodes_;
  circle_bound bound_;
  int maximum_points_;
  long level_, point_count_;
  double weight_, area_;
  bool modified_, initialized_bound_;
};

} // end namespace s2omp

#endif /* TREE_UNION_H_ */
