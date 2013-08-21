// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains a variant on the Pixel class.  The goal here is to
// use the hierarchical nature of the Pixel class to form the basis for a
// spatial quad tree structure.  Hence, a given TreePixel object will have a
// number of points associated with it, where they have been stored in such a
// way that pair finding and K nearest neighbor searches will run in ln(N) time.

#ifndef TREE_PIXEL_H_
#define TREE_PIXEL_H_

#include <map>
#include <queue>

#include "core.h"
#include "pixel.h"
#include "point.h"

namespace s2omp {
class annulus_bound; // class definition in stomp_angular_bin.h
class tree_pixel;
class tree_neighbor;
class nearest_neighbor_pixel;
class nearest_neighbor_point;

typedef struct pair_weight {
  pair_weight() : n_pairs(0), total_weight(0.0) {}
  long n_pairs;
  double total_weight;
};

typedef std::vector<tree_pixel> tree_vector;
typedef tree_vector::const_iterator tree_iterator;
typedef std::pair<tree_iterator, tree_iterator> tree_pair;
typedef std::vector<tree_pixel*> tree_ptr_vector;
typedef tree_ptr_vector::const_iterator tree_ptr_iterator;

typedef std::pair<double, tree_pixel*> distance_pixel_pair;
typedef std::priority_queue<distance_pixel_pair,
    std::vector<distance_pixel_pair>, nearest_neighbor_pixel> pixel_queue;

typedef std::pair<double, point*> distance_point_pair;
typedef std::priority_queue<distance_point_pair,
    std::vector<distance_point_pair>, nearest_neighbor_point> point_queue;

class tree_pixel: public pixel {
  // Our second variation on the pixel.  Like scalar_pixel, the idea here is
  // to use the pixel as a scaffold for sampling a field over an area.
  // Instead of storing a density, however, tree_pixel stores a vector of
  // points and the weight stored in the pixel is the sum of the weights of
  // the points.  Finally, the tree_pixel contains pointers to its
  // child-pixels.  When a point is added to the pixel, it checks the number
  // of points against the total allowed for the pixel (specified on
  // construction).  If the pixel is at capacity, it passes the point along
  // to the sub-pixels, generating a tree structure which can be traversed
  // later on for operations like pair-counting.
public:
  friend class nearest_neighbor_pixel;
  tree_pixel(uint64 id);
  tree_pixel(uint64 id, uint max_points);
  virtual ~tree_pixel();

  static tree_pixel* from_point(const point& p, int level, uint max_points);
  static tree_pixel* from_pixel(const pixel& pix, uint max_points);

  // Add a given point on the sphere to either this pixel (if the capacity for
  // this pixel hasn't been reached) or one of the sub-pixels.  Return true
  // if the point was successfully added (i.e. the point was contained in the
  // bounds of the current pixel); false, otherwise.
  bool add_point(point* p);
  bool add_point(const point& p);

  // The motivation for building a tree structure like the one in this class is
  // to do fast pair finding by recursing down the tree structure.  The work
  // of this recursion is done with _find_pairs_recursion, but the preferred
  // interfaces are via find_pairs and find_weighted_pairs, which return the
  // number of points within the input angular bin or the total weights of those
  // points.
  long find_pairs(const annulus_bound& bound) const;
  double find_weighted_pairs(const annulus_bound& bound) const;
  void _find_pairs_recursion(const annulus_bound& bound, pair_weight* pairs) const;

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
  long find_k_nearest_neighbors(
      const point& p, uint n_neighbors, point_vector* neighbors) const;

  // The special case where we're only interested in the nearest point.
  long find_nearest_neighbor(const point& p, point* neighbor) const;

  // In some cases, we're only interested in the distance to the kth nearest
  // neighbor.  The return value will be the angular distance in degrees.
  double k_nearest_neighbor_distance(
      const point& p, uint n_neighbors, long& nodes_visited) const;

  // Or in the distance to the nearest neighbor.
  double nearest_neighbor_distance(
      const point& p, long& nodes_visited) const;

  // Alternatively, we could be less interested in the nearest neighbor and
  // more interested in finding a direct match to our input point.  The
  // difference is subtle, but whereas nearest_neighbor will always return a
  // point from our tree, closest_match has an angular threshold, beyond which
  // we're not interested in the nearest neighbor because it's not a match to
  // our input point.  The returned boolean indicates whether the returned
  // point is within the specified radius and the maximum distance is
  // in degrees.
  bool closest_match(
      const point& p, double max_angular_distance, point* match) const;

  // Return the number of points contained in the current pixel and all
  // sub-pixels.
  inline long n_points() const {
    return point_count_;
  }
  inline double weight() const {
    return weight_;
  }

  // A variation on the above method, returns the number of points associated
  // with the current pixel that are also contained in the input pixel.
  long n_points(const pixel& pix) const;
  double weight(const pixel& pix) const;

  // The downside of the tree_pixel is that it doesn't really encode geometry
  // in the same way that pixels and scalar_pixels do.  This makes it hard to
  // do things like split tree_unions (defined below) into roughly equal areas
  // like we can do with pixel_unions and scalar_unions.  Coverage attempts to
  // do this based on the number of sub-nodes with data in them.  The first
  // version works on the pixel itself.  The second does the same calculation
  // for another pixel, based on the data in the current pixel.  Like the
  // unmasked fraction measures for pixels and scalar_pixels, the return values
  // cover the range [0,1].  However, the accuracy of the measure is going to
  // be a function of how many points are in the pixel (and sub-pixels) and how
  // localized they are.
  double covering_fraction() const;
  double covering_fraction(const pixel& pix) const;

  // If we want to extract a copy of all of the points that have been added
  // to this pixel, this method allows for that.
  void points(point_vector* p) const;

  // And an associated method that will extract a copy of the points associated
  // with an input pixel.
  void points(const pixel& pix, point_vector* p) const;

  // Recurse through the nodes below this one to return the number of nodes in
  // the tree.
  long n_nodes() const;

  inline uint pixel_capacity() const {
    return maximum_points_;
  }

  // Occasionally, it can be useful for outside code to be able to traverse
  // the tree structure contained in the pixel and sub-nodes.  These hooks allow
  // for access to the pointers to the sub-nodes and any point data directly.
  point_ptr_iterator points_begin() const {
    return points_.begin();
  }
  point_ptr_iterator points_end() const {
    return points_.end();
  }
  tree_ptr_iterator nodes_begin() const {
    return subnodes_.begin();
  }
  tree_ptr_iterator nodes_end() const {
    return subnodes_.end();
  }

  // And a pair of methods for indicating if the node contains points or
  // sub-nodes.
  inline bool has_nodes() const {
    return !subnodes_.empty();
  }
  inline bool has_points() const {
    return !points_.empty();
  }

  // inherited API from bound_interface
  virtual bool is_empty() const {
    return point_count_ == 0;
  }
  virtual long size() const {
    return point_count_;
  }
  virtual void clear();
  virtual double area() const {
    return exact_area() * covering_fraction();
  }

  virtual double contained_area(const pixel& pix) const {
    return exact_area() * covering_fraction(pix);
  }

private:
  tree_pixel();

  void initialize_node(uint max_points);
  bool initialize_subnodes();
  void direct_pair_count(const annulus_bound& bound, pair_weight* pairs) const;

  point_ptr_vector points_;
  tree_ptr_vector subnodes_;
  uint maximum_points_;
  long point_count_;
  double weight_;
};

class nearest_neighbor_pixel {
  // Convenience class for sorting nearest neighbor pixels in our queue.
public:
  int operator()(const distance_pixel_pair& x, const distance_pixel_pair& y) {
    // This has the opposite ordering since we want pixels ordered with the
    // closest at the top of the heap.
    return x.first > y.first;
  }
};

class nearest_neighbor_point {
  // Convenience class for sorting nearest neighbor points in our queue.
public:
  int operator()(const distance_point_pair& x, const distance_point_pair& y) {
    return x.first < y.first;
  }
};

// Default values for our nearest neigbhor finding.
static uint const DEFAULT_N_NEIGHBORS = 1;
static double const DEFAULT_MAX_NEIGHBOR_DISTANCE = 10.0;

class tree_neighbor {
  // In order to do the nearest neighbor finding in the TreePixel class, we
  // need a secondary class to handle storage of the nearest neighbor list.
  // The natural data structure for that list is a priority queue, which we're
  // using, but the fact that our list of points is sorted based on their
  // distance to a reference point means that we need a little extra plumbing
  // in order to pull that off.  Hence, the TreeNeighbor class.
public:
  friend class nearest_neighbor_point;
  tree_neighbor(const point& reference_point);
  tree_neighbor(const point& reference_point, uint n_neighbors);
  tree_neighbor(const point& reference_point, uint n_neighbors,
                double max_angular_distance);
  ~tree_neighbor();

  // Return a list of the nearest neighbors found so far.
  void nearest_neighbors(point_vector* p, bool save_neighbors);

  // Return a copy of the nearest neighbor found so far.
  point nearest_neighbor() const;

  // Return the number of neighbors in the list.  This should always be at most
  // the value used to instantiate the class, which is returned by calling
  // MaxNeighbors()
  inline uint n_neighbors() {
    return point_queue_.size();
  }
  inline uint max_neighbors() {
    return n_neighbors_;
  }

  // Submit a point for possible inclusion.  Return value indicates whether the
  // point was successfully included in the list (i.e., the distance between
  // the input point and the reference point was smaller than the current most
  // distant point in the list) or not.
  bool test_point(point* test_point);

  // Return the maximum distance of the current list.
  inline double max_distance() {
    return max_distance_;
  }

  // The default distance returned is in sin^2(theta) units since that's what
  // the edge detection code uses.  If we're interested in human units, this
  // provides that distance in degrees.
  double max_angular_distance();

  // For accounting purposes, it can be useful to keep track of how many nodes
  // we have visited during our traversal through the tree.
  inline long nodes_visited() {
    return n_nodes_visited_;
  }
  inline void add_node() {
    n_nodes_visited_++;
  }

private:
  point reference_point_;
  point_queue point_queue_;
  uint n_neighbors_;
  long n_nodes_visited_;
  double max_distance_;
};

} // end namespace s2omp

#endif /* TREE_PIXEL_H_ */
