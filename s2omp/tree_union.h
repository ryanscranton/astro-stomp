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
#include <algorithm>
#include "core.h"
#include "point.h"
#include "tree_pixel.h"
#include "pixelized_bound_interface.h"

namespace s2omp {
class angular_bin;           // class definition in stomp_angular_bin.h
class angular_correlation;   // class definition in stomp_angular_correlation.h
class radial_bin;
class pixel_union;                  // class definition in stomp_map.h
class tree_pixel;            // class definition in stomp_tree_pixel.h
class tree_map;

typedef std::map<const uint64, tree_pixel *> tree_map;
typedef tree_map::iterator tree_map_iterator;
typedef std::pair<tree_map_iterator, tree_map_iterator> tree_map_pair;

class tree_union : public bound_interface {
public:
	friend class nearest_neighbor_pixel;

	explicit tree_union(int level);
	tree_union(int level, uint16_t max_points);

	bool add_point(const point& p);

	uint32_t find_pairs(const annulus_bound& bound) const;
	void find_pairs(const point_vector& points, angular_bin* bin);

	double find_weighted_pairs(const annulus_bound& bound) const;
	void find_weighted_pairs(const point_vector& points, angular_bin* bin) const;

	// In addition to pair finding, we can also use the tree structure we've
	// built to do efficient nearest neighbor searches.  In the general case,
	// we'll be finding the k nearest neighbors of an input point.  The return
	// value is the number of nodes touched during the assemblage.
	//
	// NOTE: There is no duplication checking.  Hence, if the input point is a
	// copy of a point in the tree, then that point will be included in the
	// returned vector of points.
	uint16_t find_k_nearest_neighbors(const point& p, uint8_t n_neighbors,
			point_vector* neighbors) const;

	// The special case where we're only interested in the nearest matching point.
	uint16_t find_nearest_neighbor(const point& p, point* neighbor) const;

	// In some cases, we're only interested in the distance to the kth nearest
	// neighbor.  The return value will be the angular distance in degrees.
	double k_nearest_neighbor_distance(const point& p, uint8_t n_neighbors,
			uint16_t& nodes_visited) const;

	// Or in the distance to the nearest neighbor.
	double nearest_neighbor_distance(const point& p,
			uint16_t& nodes_visited) const;

	// Alternatively, we could be less interested in the nearest neighbor and
	// more interested in finding a direct match to our input point.  The
	// difference is subtle, but whereas NearestNeighbor will always return a
	// point from our tree, ClosestMatch has an angular threshold, beyond which
	// we're not interested in the nearest neighbor because it's not a match to
	// our input point.  The returned boolean indicates whether the returned
	// point is within the specified radius and the maximum distance is
	// in degrees.
	bool closest_match(const point& p, double max_angular_distance,
			point* match) const;

	// Return the number of points contained in the current pixel and all
	// sub-pixels.
	uint32_t n_points() const;
	double weight() const;

	// A variation on the above method, returns the number of points associated
	// with the current pixel that are also contained in the input pixel.
	uint32_t n_points(const pixel& pix) const;
	double weight(const pixel& pix) const;

	int level() const;
	uint16_t pixel_capacity() const;

	// If we want to extract a copy of all of the points that have been added
	// to this pixel, this method allows for that.
	void points(point_vector* points) const;

	// And an associated method that will extract a copy of the points associated
	// with an input pixel.
	void points(const pixel& pix, point_vector* points) const;

  // API from pixelized_bound_interface.h
	virtual bool is_empty() const;
	virtual long size() const;
	virtual void clear() const;
	virtual void area() const;

	virtual bool contains(const point& p) const;
	virtual bool contains(const pixel& pix) const;

	virtual double contained_area(const pixel& pix) const;
	virtual bool may_intersect(const pixel& pix) const;

	virtual void covering(pixel_vector* pixels) const;
	virtual void covering(int max_pixels, pixel_vector* pixels) const;
	virtual void simple_covering(int level, pixel_vector* pixels) const;

	virtual circle_bound* get_bound() const;

private:
  void neighbor_recursion(const point& p, tree_neighbor* neighbor);
  void match_recursion(const point& p, tree_neighbor* neighbor);
  void calculate_area() const;

  tree_map tree_map_;
  uint16_t maximum_points_, nodes_;
  uint32_t point_count_;
  int level_;
  double weight_, area_;
  bool modified_;
};


} // end namespace s2omp

#endif /* TREE_UNION_H_ */
