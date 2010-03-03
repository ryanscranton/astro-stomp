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

#ifndef STOMP_TREE_PIXEL_H
#define STOMP_TREE_PIXEL_H

#include <stdint.h>
#include <vector>
#include <string>
#include <map>
#include <queue>
#include "stomp_angular_coordinate.h"
#include "stomp_pixel.h"

namespace Stomp {

class AngularBin;           // class definition in stomp_angular_bin.h
class AngularCorrelation;   // class definition in stomp_angular_correlation.h
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
typedef std::priority_queue<DistancePixelPair,
  std::vector<DistancePixelPair>, NearestNeighborPixel> PixelQueue;

typedef std::pair<double, WeightedAngularCoordinate*> DistancePointPair;
typedef std::priority_queue<DistancePointPair,
  std::vector<DistancePointPair>, NearestNeighborPoint> PointQueue;

class TreePixel : public Pixel {
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
  friend class NearestNeighborPixel;
  TreePixel();
  TreePixel(const uint16_t resolution, const uint32_t pixnum,
	    const uint16_t maximum_points=200);
  TreePixel(const uint16_t resolution, const uint32_t hpixnum,
	    const uint32_t superpixnum, const uint16_t maximum_points=200);
  TreePixel(AngularCoordinate& ang, const uint16_t resolution,
	    const uint16_t maximum_points=200);
  TreePixel(const uint32_t x, const uint32_t y, const uint16_t resolution,
	    const uint16_t maximum_points=200);
  virtual ~TreePixel();

  // When the pixel has reached its carrying capacity, we want to split its
  // contents to the sub-pixels.  In this case, we create the map object that
  // contains the sub-pixels and move each of the points contained in the
  // current pixel to the correct sub-pixel.  The total weight and point count
  // for this pixel remains the same.
  bool _InitializeSubPixels();

  // The primary purpose of this class is to enable fast pair-finding for a
  // set of angular locations.  These methods implement that functionality
  // with a couple different modes of operation.  FindPairs returns an integer
  // counting of the AngularCoordinates within the specified radius or
  // annulus.  FindWeightedPairs does the same, but the value returned is
  // the sum of the weights for the objects satisfying the angular bounds.  Note
  // that the argument in this case is still an AngularCoordinate, so any
  // weight associated with that point is ignored.  The AngularCorrelation
  // versions put the number of pairs in the Counter and Weight values for each
  // angular bin.
  uint32_t DirectPairCount(AngularCoordinate& ang, AngularBin& theta,
			   int16_t region = -1);
  uint32_t FindPairs(AngularCoordinate& ang, AngularBin& theta,
		     int16_t region = -1);
  uint32_t FindPairs(AngularCoordinate& ang,
		     double theta_min, double theta_max);
  uint32_t FindPairs(AngularCoordinate& ang, double theta_max);
  double DirectWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
			     int16_t region = -1);
  double FindWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
			   int16_t region = -1);
  double FindWeightedPairs(AngularCoordinate& ang,
			   double theta_min, double theta_max);
  double FindWeightedPairs(AngularCoordinate& ang, double theta_max);

  // And for the case where we want to scale things by a weight associated with
  // each angular point explicitly.
  double DirectWeightedPairs(WeightedAngularCoordinate& w_ang,
			     AngularBin& theta, int16_t region = -1);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   AngularBin& theta, int16_t region = -1);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   double theta_min, double theta_max);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   double theta_max);

  // In these cases, the sum of the pairs are put into the Counter field for
  // the corresponding angular bin and the sum of the products of the weights
  // are put into the Weight field, if applicable.
  void FindPairs(AngularVector& ang, AngularBin& theta,
		 int16_t region = -1);
  void FindPairs(AngularVector& ang, AngularCorrelation& wtheta,
		 int16_t region = -1);
  void FindWeightedPairs(AngularVector& ang, AngularBin& theta,
			 int16_t region = -1);
  void FindWeightedPairs(AngularVector& ang, AngularCorrelation& wtheta,
			 int16_t region = -1);
  void FindWeightedPairs(WAngularVector& w_ang, AngularBin& theta,
			 int16_t region = -1);
  void FindWeightedPairs(WAngularVector& w_ang,
			 AngularCorrelation& wtheta, int16_t region = -1);

  // Since the WeightedAngularCoordinates that are fed into our tree also
  // have an arbitrary number of named Fields associated with them, we need
  // to be able to access those values as well in our pair counting.
  double DirectWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
			     const std::string& field_name,
			     int16_t region = -1);
  double FindWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
			   const std::string& field_name, int16_t region = -1);
  double FindWeightedPairs(AngularCoordinate& ang,
			   double theta_min, double theta_max,
			   const std::string& field_name);
  double FindWeightedPairs(AngularCoordinate& ang, double theta_max,
			   const std::string& field_name);
  void FindWeightedPairs(AngularVector& ang, AngularBin& theta,
			 const std::string& field_name, int16_t region = -1);
  void FindWeightedPairs(AngularVector& ang, AngularCorrelation& wtheta,
			 const std::string& field_name, int16_t region = -1);

  // If we have a WeightedAngularCoordinate as the input, then we need to
  // account for the case where you want to use the weight associated with
  // the coordinate as well as the case where you want to use a field from
  // the input coordinate.  First the Weight vs. Field case.
  double DirectWeightedPairs(WeightedAngularCoordinate& w_ang,
			     AngularBin& theta, const std::string& field_name,
			     int16_t region = -1);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang, AngularBin& theta,
			   const std::string& field_name, int16_t region = -1);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   double theta_min, double theta_max,
			   const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang, double theta_max,
			   const std::string& field_name);
  void FindWeightedPairs(WAngularVector& w_ang, AngularBin& theta,
			 const std::string& field_name, int16_t region = -1);
  void FindWeightedPairs(WAngularVector& w_ang, AngularCorrelation& wtheta,
			 const std::string& field_name, int16_t region = -1);

  // And finally, the Field vs. Field case.
  double DirectWeightedPairs(WeightedAngularCoordinate& w_ang,
			     const std::string& ang_field_name,
			     AngularBin& theta, const std::string& field_name,
			     int16_t region = -1);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   const std::string& ang_field_name, AngularBin& theta,
			   const std::string& field_name, int16_t region = -1);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   const std::string& ang_field_name,
			   double theta_min, double theta_max,
			   const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   const std::string& ang_field_name, double theta_max,
			   const std::string& field_name);
  void FindWeightedPairs(WAngularVector& w_ang,
			 const std::string& ang_field_name,
			 AngularBin& theta, const std::string& field_name,
			 int16_t region = -1);
  void FindWeightedPairs(WAngularVector& w_ang,
			 const std::string& ang_field_name,
			 AngularCorrelation& wtheta,
			 const std::string& field_name, int16_t region = -1);

  // In addition to pair finding, we can also use the tree structure we've
  // built to do efficient nearest neighbor searches.  In the general case,
  // we'll be finding the k nearest neighbors of an input point.  The return
  // value is the number of nodes touched during the assemblage.
  //
  // NOTE: There is no duplication checking.  Hence, if the input point is a
  // copy of a point in the tree, then that point will be included in the
  // returned vector of points.
  uint16_t FindKNearestNeighbors(AngularCoordinate& ang, uint8_t n_neighbors,
				 WAngularVector& neighbors_ang);

  // The special case where we're only interested in the nearest matching point.
  uint16_t FindNearestNeighbor(AngularCoordinate& ang,
			       WeightedAngularCoordinate& neighbor_ang);

  // In some cases, we're only interested in the distance to the kth nearest
  // neighbor.  The return value will be the angular distance in degrees.
  double KNearestNeighborDistance(AngularCoordinate& ang, uint8_t n_neighbors,
				  uint16_t& nodes_visited);

  // Or in the distance to the nearest neighbor.
  double NearestNeighborDistance(AngularCoordinate& ang,
				 uint16_t& nodes_visited);

  // For the recursion necessary to do the neighbor finding, we use this
  // internal method.
  void _NeighborRecursion(AngularCoordinate& ang, TreeNeighbor& neighbor);

  // For the pair finding, we end up checking the X-Y-Z corners of the pixels
  // a lot, so we store those values internally and use an internal method for
  // finding the intersection of the pixel with the input annulus;
  virtual int8_t IntersectsAnnulus(AngularCoordinate& ang, AngularBin& theta);

  // We also have an internal version of the code for finding the edge
  // distances.  The return boolean tells us whether the distance returned
  // is to a pixel edge (true) or a pixel corner (false).
  virtual bool EdgeDistances(AngularCoordinate& ang, double& min_edge_distance,
			     double& max_edge_distance);

  // And a method to set these values up internally.
  void InitializeCorners();

  // Add a given point on the sphere to either this pixel (if the capacity for
  // this pixel hasn't been reached) or one of the sub-pixels.  Return true
  // if the point was successfully added (i.e. the point was contained in the
  // bounds of the current pixel); false, otherwise.
  bool AddPoint(WeightedAngularCoordinate* ang);

  // The default method for adding WeightedAngularCoordinates to the pixel
  // takes a pointer to the object.  This means that the pixel now owns that
  // object and it shouldn't be deleted from the heap except by the pixel.
  // For cases where we want to retain a copy of the point outside of the
  // pixel, we provide a second method which takes a reference to the object
  // and creates and stores an internal copy.  The input object can thus be
  // modified or deleted without affecting the tree.
  bool AddPoint(WeightedAngularCoordinate& w_ang);

  // Complimentary method for specifying a weight separately when adding a
  // point to the pixel.
  bool AddPoint(AngularCoordinate& ang, double object_weight = 1.0);

  // Return the number of points contained in the current pixel and all
  // sub-pixels.
  uint32_t NPoints();

  // A variation on the above method, returns the number of points associated
  // with the current pixel that are also contained in the input pixel.
  uint32_t NPoints(Pixel& pix);

  // Likewise, we can provide a similar method for returning the weight
  // associated with an input pixel.
  double PixelWeight(Pixel& pix);

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
  void Points(WAngularVector& w_ang);

  // And an associated method that will extract a copy of the points associated
  // with an input pixel.
  void Points(WAngularVector& w_ang, Pixel& pix);

  // Recurse through the nodes below this one to return the number of nodes in
  // the tree.
  uint16_t Nodes();
  void _AddSubNodes(uint16_t& n_nodes);

  // Modify the weight of the pixel.  Generally this is only called when adding
  // a point to the pixel.  Calling it directly will result in a pixel weight
  // which is no longer the sum of the contained points' weights.
  void AddToWeight(double weight);

  // Since our WeightedAngularCoordinate objects have an arbitrary number
  // of Fields associated with them, we store that information as well when
  // we're building our tree structure.  These methods allow for access to
  // the aggregate values for a given Field.
  double FieldTotal(const std::string& field_name);
  double FieldTotal(const std::string& field_name, Pixel& pix);
  void AddToField(const std::string& field_name, double weight);
  uint16_t NField();
  bool HasFields();
  void FieldNames(std::vector<std::string>& field_names);

  // Modify and return the point capacity for the pixel, respectively.
  void SetPixelCapacity(uint16_t maximum_points);
  uint16_t PixelCapacity();

  // Since we're storing pointers to the WeightedAngularCoordinates, we need
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
  WAngularPtrVector ang_;
  FieldDict field_total_;
  uint16_t maximum_points_;
  uint32_t point_count_;
  bool initialized_subpixels_;
  double unit_sphere_x_, unit_sphere_y_, unit_sphere_z_;
  double unit_sphere_x_ul_, unit_sphere_y_ul_, unit_sphere_z_ul_;
  double unit_sphere_x_ll_, unit_sphere_y_ll_, unit_sphere_z_ll_;
  double unit_sphere_x_ur_, unit_sphere_y_ur_, unit_sphere_z_ur_;
  double unit_sphere_x_lr_, unit_sphere_y_lr_, unit_sphere_z_lr_;
  TreePtrVector subpix_;
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
  TreeNeighbor(AngularCoordinate& reference_ang,
	       uint8_t n_neighbors = 1);
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

} // end namespace Stomp

#endif
