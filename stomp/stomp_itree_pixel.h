// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains a variant on the TreePixel class.  Again, we are
// using the hierarchical pixelization as the scaffolding for a spatial tree
// structure.  However, in this case, the tree is populated by
// IndexedAngularCoordinate objects rather than their weighted cousins.  For
// integer data types, we're not going to be interested in aggregate statistics
// like with the floating point data, but rather returning vectors of objects
// or indices.


#ifndef STOMP_ITREE_PIXEL_H
#define STOMP_ITREE_PIXEL_H

#include <stdint.h>
#include <vector>
#include <string>
#include <map>
#include <queue>
#include "stomp_angular_coordinate.h"
#include "stomp_pixel.h"

namespace Stomp {

class AngularBin;           // class definition in stomp_angular_bin.h
class IndexedTreePixel;
class IndexedTreeNeighbor;
class NearestNeighborIndexedPixel;
class NearestNeighborIndexedPoint;

typedef std::vector<IndexedTreePixel> ITreeVector;
typedef ITreeVector::iterator ITreeIterator;
typedef std::pair<ITreeIterator, ITreeIterator> ITreePair;
typedef std::vector<IndexedTreePixel *> ITreePtrVector;
typedef ITreePtrVector::iterator ITreePtrIterator;

typedef std::pair<double, IndexedTreePixel*> DistanceIPixelPair;
typedef std::priority_queue<DistanceIPixelPair,
  std::vector<DistanceIPixelPair>, NearestNeighborIndexedPixel> IPixelQueue;

typedef std::pair<double, IndexedAngularCoordinate*> DistanceIPointPair;
typedef std::priority_queue<DistanceIPointPair,
  std::vector<DistanceIPointPair>, NearestNeighborIndexedPoint> IPointQueue;

typedef std::vector<uint32_t> IndexVector;
typedef IndexVector::iterator IndexIterator;

class IndexedTreePixel : public Pixel {
  // Similiar in spirit to the TreePixel, except that we're storing
  // IndexedAngularCoordinates.  The Weight() field we inherit from the Pixel
  // class will be meaningless in this case, however, as we'll be mostly
  // interested in returning copies of the stored objects or their indices.
 public:
  friend class NearestNeighborIndexedPixel;
  IndexedTreePixel();
  IndexedTreePixel(const uint32_t resolution, const uint32_t pixnum,
		   const uint16_t maximum_points=200);
  IndexedTreePixel(AngularCoordinate& ang, const uint32_t resolution,
		   const uint16_t maximum_points=200);
  IndexedTreePixel(const uint32_t x, const uint32_t y,
		   const uint32_t resolution,
		   const uint16_t maximum_points=200);
  virtual ~IndexedTreePixel();

  // When the pixel has reached its carrying capacity, we want to split its
  // contents to the sub-pixels.  In this case, we create the map object that
  // contains the sub-pixels and move each of the points contained in the
  // current pixel to the correct sub-pixel.
  bool _InitializeSubPixels();

  // The primary purpose of this class is to enable fast pair-finding for a
  // set of angular locations.  These methods implement that functionality
  // with a couple different modes of operation.  FindPairs returns either a
  // vector of copies of the tree objects within the angular bounds or a
  // vector of the index values for those objects.  Unlike the TreePixel
  // variants, the AngularBin objects used in some of the methods do not store
  // any pair-finding results.
  void FindPairs(AngularCoordinate& ang, AngularBin& theta,
		 IAngularVector& i_angVec);
  void FindPairs(AngularCoordinate& ang, AngularBin& theta,
		 IndexVector& pair_indices);
  void FindPairs(AngularCoordinate& ang,
		 double theta_min, double theta_max,
		 IAngularVector& i_angVec);
  void FindPairs(AngularCoordinate& ang,
		 double theta_min, double theta_max,
		 IndexVector& pair_indices);
  void FindPairs(AngularCoordinate& ang, double theta_max,
		 IAngularVector& pair_indices);
  void FindPairs(AngularCoordinate& ang, double theta_max,
		 IndexVector& pair_indices);

  // The above methods should be called for pair-finding on the tree objects.
  // The following two internal methods handle the recursion within the
  // tree and should not be called directly.
  void _FindPairs(AngularCoordinate& ang, AngularBin& theta,
		  IAngularPtrVector& i_ang);
  void _PointPtrs(IAngularPtrVector& i_ang);

  // In addition to pair finding, we can also use the tree structure we've
  // built to do efficient nearest neighbor searches.  In the general case,
  // we'll be finding the k nearest neighbors of an input point.  The return
  // value is the number of nodes touched during the assemblage.
  //
  // NOTE: There is no duplication checking.  Hence, if the input point is a
  // copy of a point in the tree, then that point will be included in the
  // returned vector of points.
  uint16_t FindKNearestNeighbors(AngularCoordinate& ang, uint8_t n_neighbors,
				 IAngularVector& neighbors_ang);

  // The special case where we're only interested in the nearest matching point.
  uint16_t FindNearestNeighbor(AngularCoordinate& ang,
			       IndexedAngularCoordinate& neighbor_ang);

  // In some cases, we're only interested in the distance to the kth nearest
  // neighbor.  The return value will be the angular distance in degrees.
  double KNearestNeighborDistance(AngularCoordinate& ang, uint8_t n_neighbors,
				  uint16_t& nodes_visited);

  // Or in the distance to the nearest neighbor.
  double NearestNeighborDistance(AngularCoordinate& ang,
				 uint16_t& nodes_visited);

  // Alternatively, we could be less interested in the nearest neighbor and
  // more interested in finding a direct match to our input point.  The
  // difference is subtle, but whereas NearestNeighbor will always return a
  // point from our tree, ClosestMatch has an angular threshold, beyond which
  // we're not interested in the nearest neighbor because it's not a match to
  // our input point.  The returned boolean indicates whether the returned
  // point is within the specified radius and the maximum distance is
  // in degrees.
  bool ClosestMatch(AngularCoordinate& ang, double max_distance,
		    IndexedAngularCoordinate& match_ang);

  // For the recursion necessary to do the neighbor finding, we use this
  // internal method.
  void _NeighborRecursion(AngularCoordinate& ang,
			  IndexedTreeNeighbor& neighbor);

  // And a method to set these values up internally.
  void InitializeCorners();

  // Add a given point on the sphere to either this pixel (if the capacity for
  // this pixel hasn't been reached) or one of the sub-pixels.  Return true
  // if the point was successfully added (i.e. the point was contained in the
  // bounds of the current pixel); false, otherwise.
  bool AddPoint(IndexedAngularCoordinate* ang);

  // The default method for adding IndexedAngularCoordinates to the pixel
  // takes a pointer to the object.  This means that the pixel now owns that
  // object and it shouldn't be deleted from the heap except by the pixel.
  // For cases where we want to retain a copy of the point outside of the
  // pixel, we provide a second method which takes a reference to the object
  // and creates and stores an internal copy.  The input object can thus be
  // modified or deleted without affecting the tree.
  bool AddPoint(IndexedAngularCoordinate& w_ang);

  // Complimentary method for specifying a weight separately when adding a
  // point to the pixel.
  bool AddPoint(AngularCoordinate& ang, uint32_t index);

  // Return the number of points contained in the current pixel and all
  // sub-pixels.
  uint32_t NPoints();

  // A variation on the above method, returns the number of points associated
  // with the current pixel that are also contained in the input pixel.
  uint32_t NPoints(Pixel& pix);

  // Likewise, we can provide a similar method for returning the indices
  // associated with an input pixel.
  void Indices(Pixel& pix, IndexVector& indices);

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
  double Coverage();
  double Coverage(Pixel& pix);

  // If we want to extract a copy of all of the points that have been added
  // to this pixel, this method allows for that.
  void Points(IAngularVector& i_ang);

  // And an associated method that will extract a copy of the points associated
  // with an input pixel.
  void Points(IAngularVector& i_ang, Pixel& pix);

  // Recurse through the nodes below this one to return the number of nodes in
  // the tree.
  uint16_t Nodes();
  void _AddSubNodes(uint16_t& n_nodes);

  // Modify and return the point capacity for the pixel, respectively.
  void SetPixelCapacity(uint16_t maximum_points);
  uint16_t PixelCapacity();

  // Occasionally, it can be useful for outside code to be able to traverse
  // the tree structure contained in the pixel and sub-nodes.  These hooks allow
  // for access to the pointers to the sub-nodes and any point data directly.
  IAngularPtrIterator PointsBegin();
  IAngularPtrIterator PointsEnd();
  ITreePtrIterator NodesBegin();
  ITreePtrIterator NodesEnd();

  // And a pair of methods for indicating if the node contains points or
  // sub-nodes.
  bool HasPoints();
  bool HasNodes();

  // Since we're storing pointers to the IndexedAngularCoordinates, we need
  // to explicitly delete them to clear all of the memory associated with the
  // pixel.
  void Clear();

  // Since we've got this data stored locally in variables, we can use faster
  // accessors than the standard Pixel methods.
  virtual double UnitSphereX();
  virtual double UnitSphereY();
  virtual double UnitSphereZ();

  virtual double UnitSphereX_UL();
  virtual double UnitSphereY_UL();
  virtual double UnitSphereZ_UL();

  virtual double UnitSphereX_UR();
  virtual double UnitSphereY_UR();
  virtual double UnitSphereZ_UR();

  virtual double UnitSphereX_LL();
  virtual double UnitSphereY_LL();
  virtual double UnitSphereZ_LL();

  virtual double UnitSphereX_LR();
  virtual double UnitSphereY_LR();
  virtual double UnitSphereZ_LR();

  // We can also speed up some of the other pixel boundary checking routines
  // by using internal versions.
  virtual void WithinAnnulus(AngularBin& theta, PixelVector& pix,
			     bool check_full_pixel);

 private:
  IAngularPtrVector ang_;
  uint16_t maximum_points_;
  uint32_t point_count_;
  bool initialized_subpixels_;
  double unit_sphere_x_, unit_sphere_y_, unit_sphere_z_;
  double unit_sphere_x_ul_, unit_sphere_y_ul_, unit_sphere_z_ul_;
  double unit_sphere_x_ll_, unit_sphere_y_ll_, unit_sphere_z_ll_;
  double unit_sphere_x_ur_, unit_sphere_y_ur_, unit_sphere_z_ur_;
  double unit_sphere_x_lr_, unit_sphere_y_lr_, unit_sphere_z_lr_;
  ITreePtrVector subpix_;
};

class NearestNeighborIndexedPixel {
  // Convenience class for sorting nearest neighbor pixels in our queue.
 public:
  int operator()(const DistanceIPixelPair& x, const DistanceIPixelPair& y) {
    // This has the opposite ordering since we want pixels ordered with the
    // closest at the top of the heap.
    return x.first > y.first;
  }
};

class NearestNeighborIndexedPoint {
  // Convenience class for sorting nearest neighbor points in our queue.
 public:
  int operator()(const DistanceIPointPair& x, const DistanceIPointPair& y) {
    return x.first < y.first;
  }
};

class IndexedTreeNeighbor {
  // In order to do the nearest neighbor finding in the IndexedTreePixel
  // class, we need a secondary class to handle storage of the nearest neighbor
  // list.  The natural data structure for that list is a priority queue, which
  // we're using, but the fact that our list of points is sorted based on their
  // distance to a reference point means that we need a little extra plumbing
  // in order to pull that off.  Hence, the TreeNeighbor class.
 public:
  friend class NearestNeighborIndexedPoint;
  IndexedTreeNeighbor(AngularCoordinate& reference_ang,
	       uint8_t n_neighbors = 1);
  IndexedTreeNeighbor(AngularCoordinate& reference_ang,
	       uint8_t n_neighbors, double max_distance);
  ~IndexedTreeNeighbor();

  // Return a list of the nearest neighbors found so far.
  void NearestNeighbors(IAngularVector& i_ang, bool save_neighbors = true);

  // Return the number of neighbors in the list.  This should always be at most
  // the value used to instantiate the class, which is returned by calling
  // MaxNeighbors()
  uint8_t Neighbors();
  uint8_t MaxNeighbors();

  // Submit a point for possible inclusion.  Return value indicates whether the
  // point was successfully included in the list (i.e., the distance between
  // the input point and the reference point was smaller than the current most
  // distant point in the list) or not.
  bool TestPoint(IndexedAngularCoordinate* test_ang);

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
  IPointQueue ang_queue_;
  uint8_t n_neighbors_;
  uint16_t n_nodes_visited_;
  double max_distance_;
};

} // end namespace Stomp

#endif
