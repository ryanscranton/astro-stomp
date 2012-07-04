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

namespace s2omp {
class AngularBin; // class definition in stomp_angular_bin.h
class AngularCorrelation; // class definition in stomp_angular_correlation.h
class RadialBin;
class TreePixel;
class TreeNeighbor;
class NearestNeighborPixel;
class NearestNeighborPoint;

typedef std::vector<TreePixel> TreeVector;
typedef TreeVector::iterator TreeIterator;
typedef std::pair<TreeIterator, TreeIterator> TreePair;
typedef std::vector<TreePixel *> TreePtrVector;
typedef TreePtrVector::iterator TreePtrIterator;

typedef std::pair<double, TreePixel*> DistancePixelPair;
typedef std::priority_queue<DistancePixelPair, std::vector<DistancePixelPair>,
		NearestNeighborPixel> PixelQueue;

typedef std::pair<double, WeightedAngularCoordinate*> DistancePointPair;
typedef std::priority_queue<DistancePointPair, std::vector<DistancePointPair>,
		NearestNeighborPoint> PointQueue;

class tree_pixel: public pixel {
	// Our second variation on the Pixel.  Like ScalarPixel, the idea
	// here is to use the Pixel as a scaffold for sampling a field over an
	// area.  Instead of storing a density, however, TreePixel stores a
	// vector of WeightedAngularCoordinates and the weight stored in the pixel
	// is the sum of the weights of the WeightedAngularCoordinates.  Finally,
	// the TreePixel contains pointers to its sub-pixels.  When a point is
	// added to the pixel, it checks the number of points against the total
	// allowed for the pixel (specified on construction).  If the pixel is at
	// capacity, it passes the point along to the sub-pixels, generating a tree
	// structure which can be traversed later on for operations like
	// pair-counting.
public:
	friend class nearest_neighbor_pixel;
	explicit tree_pixel(uint64 id);
	tree_pixel(uint64 id, uint16_t max_points);
	virtual ~tree_pixel();

	static tree_pixel from_point(const point& p, int level,
			uint16_t maximum_points);
	static tree_pixel from_pixel(const pixel& pix, uint16_t maximum_points);

	// Add a given point on the sphere to either this pixel (if the capacity for
	// this pixel hasn't been reached) or one of the sub-pixels.  Return true
	// if the point was successfully added (i.e. the point was contained in the
	// bounds of the current pixel); false, otherwise.
	bool add_point(const point& p);

	uint32_t find_pairs(const annulus_bound& bound) const;
	double find_weighted_pairs(const annulus_bound& bound) const;

	// In addition to pair finding, we can also use the tree structure we've
	// built to do efficient nearest neighbor searches.  In the general case,
	// we'll be finding the k nearest neighbors of an input point.  The return
	// value is the number of nodes touched during the assemblage.
	//
	// NOTE: There is no duplication checking.  Hence, if the input point is a
	// copy of a point in the tree, then that point will be included in the
	// returned vector of points.
	uint16_t find_k_nearest_neighbors(const point& p, uint8_t n_neighbors,
			point_vector& neighbors) const;

	// The special case where we're only interested in the nearest matching point.
	uint16_t find_nearest_neighbor(const point& p, point& neighbor) const;

	// In some cases, we're only interested in the distance to the kth nearest
	// neighbor.  The return value will be the angular distance in degrees.
	double k_nearest_neighbor_distance(const point& p, uint8_t n_neighbors,
			uint16_t& nodes_visited) const;

	// Or in the distance to the nearest neighbor.
	double
			nearest_neighbor_distance(const point& p, uint16_t& nodes_visited) const;

	// Alternatively, we could be less interested in the nearest neighbor and
	// more interested in finding a direct match to our input point.  The
	// difference is subtle, but whereas NearestNeighbor will always return a
	// point from our tree, ClosestMatch has an angular threshold, beyond which
	// we're not interested in the nearest neighbor because it's not a match to
	// our input point.  The returned boolean indicates whether the returned
	// point is within the specified radius and the maximum distance is
	// in degrees.
	bool
			closest_match(const point& p, double max_angular_distance, point& match) const;

	// Return the number of points contained in the current pixel and all
	// sub-pixels.
	uint32_t n_points() const;
	double weight() const;

	// A variation on the above method, returns the number of points associated
	// with the current pixel that are also contained in the input pixel.
	uint32_t n_points(const pixel& pix) const;
	double weight(const pixel& pix) const;

	// The downside of the TreePixel is that it doesn't really encode geometry
	// in the same way that Pixels and ScalarPixels do.  This makes it hard to
	// do things like split TreeMaps (defined below) into roughly equal areas
	// like we can do with Maps and ScalarMaps.  Coverage attempts to do this
	// based on the number of sub-nodes with data in them.  The first version
	// works on the pixel itself.  The second does the same calculation for
	// another pixel, based on the data in the current pixel.  Like the unmasked
	// fraction measures for Pixels and ScalarPixels, the return values cover
	// the range [0,1].  However, the accuracy of the measure is going to be a
	// function of how many points are in the pixel (and sub-pixels) and how
	// localized they are.
	double covering_fraction() const;
	double covering_fraction(const pixel& pix) const;

	// If we want to extract a copy of all of the points that have been added
	// to this pixel, this method allows for that.
	void points(point_vector& points) const;

	// And an associated method that will extract a copy of the points associated
	// with an input pixel.
	void points(const pixel& pix, point_vector& points) const;

	// Recurse through the nodes below this one to return the number of nodes in
	// the tree.
	uint16_t n_nodes() const;

	uint16_t pixel_capacity() const;

	// Occasionally, it can be useful for outside code to be able to traverse
	// the tree structure contained in the pixel and sub-nodes.  These hooks allow
	// for access to the pointers to the sub-nodes and any point data directly.
	point_ptr_iterator points_begin();
	point_ptr_iterator points_end();
	tree_ptr_iterator nodes_begin();
	tree_ptr_iterator nodes_end();

	// And a pair of methods for indicating if the node contains points or
	// sub-nodes.
	bool has_points() const;
	bool has_nodes() const;

	// Since we're storing pointers to the WeightedAngularCoordinates, we need
	// to explicitly delete them to clear all of the memory associated with the
	// pixel.
	void clear();

private:
	tree_pixel();

	bool _initialize_children();
	void _add_children(uint16_t& n_nodes);
	uint32_t direct_pair_count(annulus_bound& bound);
	double direct_weighted_pairs(annulus_bound& bound);
	void _neighbor_recursion(point& p, tree_neighbor& neighbor);

	point_ptr_vector points_;
	uint16_t maximum_points_;
	uint32_t point_count_;
	double weight_;
	bool initialized_children_;
	tree_ptr_vector children_;
};

class NearestNeighborPixel {
	// Convenience class for sorting nearest neighbor pixels in our queue.
public:
	int operator()(const DistancePixelPair& x, const DistancePixelPair& y) {
		// This has the opposite ordering since we want pixels ordered with the
		// closest at the top of the heap.
		return x.first > y.first;
	}
};

class NearestNeighborPoint {
	// Convenience class for sorting nearest neighbor points in our queue.
public:
	int operator()(const DistancePointPair& x, const DistancePointPair& y) {
		return x.first < y.first;
	}
};

class TreeNeighbor {
	// In order to do the nearest neighbor finding in the TreePixel class, we
	// need a secondary class to handle storage of the nearest neighbor list.
	// The natural data structure for that list is a priority queue, which we're
	// using, but the fact that our list of points is sorted based on their
	// distance to a reference point means that we need a little extra plumbing
	// in order to pull that off.  Hence, the TreeNeighbor class.
public:
	friend class NearestNeighborPoint;
	TreeNeighbor(AngularCoordinate& reference_ang, uint8_t n_neighbors = 1);
	TreeNeighbor(AngularCoordinate& reference_ang, uint8_t n_neighbors,
			double max_distance);
	~TreeNeighbor();

	// Return a list of the nearest neighbors found so far.
	void NearestNeighbors(WAngularVector& w_ang, bool save_neighbors = true);

	// Return the number of neighbors in the list.  This should always be at most
	// the value used to instantiate the class, which is returned by calling
	// MaxNeighbors()
	uint8_t Neighbors();
	uint8_t MaxNeighbors();

	// Submit a point for possible inclusion.  Return value indicates whether the
	// point was successfully included in the list (i.e., the distance between
	// the input point and the reference point was smaller than the current most
	// distant point in the list) or not.
	bool TestPoint(WeightedAngularCoordinate* test_ang);

	// Return the maximum distance of the current list.
	double MaxDistance();

	// The default distance returned is in sin^2(theta) units since that's what
	// the edge detection code uses.  If we're interested in human units, this
	// provides that distance in degrees.
	double MaxAngularDistance();

	// For accounting purposes, it can be useful to keep track of how many nodes
	// we have visited during our traversal through the tree.
	uint16_t NodesVisited();
	void AddNode();

private:
	AngularCoordinate reference_ang_;
	PointQueue ang_queue_;
	uint8_t n_neighbors_;
	uint16_t n_nodes_visited_;
	double max_distance_;
};

} // end namespace s2omp

#endif /* TREE_PIXEL_H_ */
