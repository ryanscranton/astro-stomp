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

#ifndef STOMP_TREE_MAP_H
#define STOMP_TREE_MAP_H

#include <stdint.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <algorithm>
#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "stomp_tree_pixel.h"
#include "stomp_base_map.h"

namespace Stomp {

class AngularBin;           // class definition in stomp_angular_bin.h
class AngularCorrelation;   // class definition in stomp_angular_correlation.h
class Map;                  // class definition in stomp_map.h
class TreePixel;            // class definition in stomp_tree_pixel.h
class TreeMap;

typedef std::map<const uint32_t, TreePixel *> TreeDict;
typedef TreeDict::iterator TreeDictIterator;
typedef std::pair<TreeDictIterator, TreeDictIterator> TreeDictPair;

typedef std::vector<TreeMap> TreeMapVector;
typedef TreeMapVector::iterator TreeMapIterator;
typedef std::pair<TreeMapIterator, TreeMapIterator> TreeMapPair;

class TreeMap : public BaseMap {
  // Another variation on the Map.  Unlike the Map and DensityMap
  // classes, TreeMap objects start off with no defined geometry.  Instead,
  // as AngularCoordinate objects are placed into the TreeMap instance,
  // the map automatically generates the necessary nodes for storing the data.
  // The resulting structure can then be used to very quickly find object
  // pairs on a number of different scales.

 public:
  friend class NearestNeighborPixel;
  // Since the geometry is specified when points are added to the map, the
  // number of parameters for the constructor is potentially very small.  By
  // default the constructor will set things up to have the base resolution of
  // the tree match the resolution of the superpixels and the maximum number of
  // points per node to be 50.  The latter is somewhat arbitrary.  However, the
  // resolution for the base level of the map is important if you want to do
  // pair counting in conjunction with a DensityMap that uses sub-regions.
  // The region resolution for that DensityMap must match the resolution
  // chosen for the base level of the TreeMap.
  TreeMap(uint32_t resolution=HPixResolution, uint16_t maximum_points=50);

  // Alternatively, we can make use of the Read method (see below) to
  // populate our TreeMap during instantion, given an ascii file containing
  // angular coordinate data.
  TreeMap(const std::string& input_file,
	  uint32_t resolution=HPixResolution, uint16_t maximum_points=50,
	  AngularCoordinate::Sphere sphere = AngularCoordinate::Equatorial,
	  bool verbose = false, uint8_t theta_column = 0,
	  uint8_t phi_column = 1, int8_t weight_column = -1);
  TreeMap(const std::string& input_file, FieldColumnDict& field_columns,
	  uint32_t resolution=HPixResolution, uint16_t maximum_points=50,
	  AngularCoordinate::Sphere sphere = AngularCoordinate::Equatorial,
	  bool verbose = false, uint8_t theta_column = 0,
	  uint8_t phi_column = 1, int8_t weight_column = -1);
  ~TreeMap();

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
  uint32_t FindPairs(AngularCoordinate& ang, AngularBin& theta);
  uint32_t FindPairs(AngularCoordinate& ang,
		     double theta_min, double theta_max);
  uint32_t FindPairs(AngularCoordinate& ang, double theta_max);
  void FindPairs(AngularVector& ang, AngularBin& theta);
  void FindPairs(AngularVector& ang, AngularCorrelation& wtheta);

  double FindWeightedPairs(AngularCoordinate& ang, AngularBin& theta);
  double FindWeightedPairs(AngularCoordinate& ang,
			   double theta_min, double theta_max);
  double FindWeightedPairs(AngularCoordinate& ang, double theta_max);
  void FindWeightedPairs(AngularVector& ang, AngularBin& theta);
  void FindWeightedPairs(AngularVector& ang, AngularCorrelation& wtheta);

  // And for the case where we want to scale things by a weight associated with
  // each angular point explicitly.
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   AngularBin& theta);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   double theta_min, double theta_max);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   double theta_max);
  void FindWeightedPairs(WAngularVector& w_ang, AngularBin& theta);
  void FindWeightedPairs(WAngularVector& w_ang,
			 AngularCorrelation& wtheta);

  // And for the cases where we want to access the Field values in the tree.
  double FindWeightedPairs(AngularCoordinate& ang, AngularBin& theta,
			   const std::string& field_name);
  double FindWeightedPairs(AngularCoordinate& ang,
			   double theta_min, double theta_max,
			   const std::string& field_name);
  double FindWeightedPairs(AngularCoordinate& ang, double theta_max,
			   const std::string& field_name);
  void FindWeightedPairs(AngularVector& ang, AngularBin& theta,
			 const std::string& field_name);
  void FindWeightedPairs(AngularVector& ang, AngularCorrelation& wtheta,
			 const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang, AngularBin& theta,
			   const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   double theta_min, double theta_max,
			   const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang, double theta_max,
			   const std::string& field_name);
  void FindWeightedPairs(WAngularVector& w_ang, AngularBin& theta,
			 const std::string& field_name);
  void FindWeightedPairs(WAngularVector& w_ang, AngularCorrelation& wtheta,
			 const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   const std::string& ang_field_name, AngularBin& theta,
			   const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   const std::string& ang_field_name,
			   double theta_min, double theta_max,
			   const std::string& field_name);
  double FindWeightedPairs(WeightedAngularCoordinate& w_ang,
			   const std::string& ang_field_name, double theta_max,
			   const std::string& field_name);
  void FindWeightedPairs(WAngularVector& w_ang,
			 const std::string& ang_field_name,
			 AngularBin& theta, const std::string& field_name);
  void FindWeightedPairs(WAngularVector& w_ang,
			 const std::string& ang_field_name,
			 AngularCorrelation& wtheta,
			 const std::string& field_name);

  // And for a selected set of the above variations, we also include forms
  // which allow the regions to come into play.  Generally speaking, these
  // are the versions that are most likely to be called from the correlation
  // codes.
  void FindPairsWithRegions(AngularVector& ang, AngularBin& theta);
  void FindPairsWithRegions(AngularVector& ang, AngularCorrelation& wtheta);
  void FindWeightedPairsWithRegions(AngularVector& ang, AngularBin& theta);
  void FindWeightedPairsWithRegions(AngularVector& ang,
                                    AngularCorrelation& wtheta);
  void FindWeightedPairsWithRegions(WAngularVector& w_ang, AngularBin& theta);
  void FindWeightedPairsWithRegions(WAngularVector& w_ang,
                                    AngularCorrelation& wtheta);
  void FindWeightedPairsWithRegions(AngularVector& ang, AngularBin& theta,
                                    const std::string& field_name);
  void FindWeightedPairsWithRegions(AngularVector& ang,
                                    AngularCorrelation& wtheta,
                                    const std::string& field_name);
  void FindWeightedPairsWithRegions(WAngularVector& w_ang, AngularBin& theta,
                                    const std::string& field_name);
  void FindWeightedPairsWithRegions(WAngularVector& w_ang,
                                    AngularCorrelation& wtheta,
                                    const std::string& field_name);
  void FindWeightedPairsWithRegions(WAngularVector& w_ang,
                                    const std::string& ang_field_name,
                                    AngularBin& theta,
                                    const std::string& field_name);
  void FindWeightedPairsWithRegions(WAngularVector& w_ang,
                                    const std::string& ang_field_name,
                                    AngularCorrelation& wtheta,
                                    const std::string& field_name);

  // In addition to pair finding, we can also use the tree structure we've
  // built to do efficient nearest neighbor searches.  In the general case,
  // we'll be finding the k nearest neighbors of an input point.  The return
  // value is the number of nodes visited during assemblage.
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


  // Add a given point on the sphere to the map.
  bool AddPoint(WeightedAngularCoordinate* ang);

  // The default method for adding WeightedAngularCoordinates to the map
  // takes a pointer to the object.  This means that the map now owns that
  // object and it shouldn't be deleted from the heap except by the map.
  // For cases where we want to retain a copy of the point outside of the
  // map, we provide a second method which takes a reference to the object
  // and creates and stores an internal copy.  The input object can thus be
  // modified or deleted without affecting the map.
  bool AddPoint(WeightedAngularCoordinate& w_ang);

  // Complimentary method for specifying a weight separately when adding a
  // point to the pixel.
  bool AddPoint(AngularCoordinate& ang, double object_weight = 1.0);

  // Rather than adding points one by one, we can also take an input file and
  // add those points to the tree.  We can do this with and without also adding
  // Field values to each point from the input file.  If the weight column is
  // not specified, then the Weight is set to unity for each point.  As with
  // the ToWAngularVector methods in the WeightedAngularCoordinate class, the
  // returned boolean indicates success or failure in adding all of the points.
  bool Read(const std::string& input_file,
	    AngularCoordinate::Sphere sphere = AngularCoordinate::Equatorial,
	    bool verbose = false, uint8_t theta_column = 0,
	    uint8_t phi_column = 1, int8_t weight_column = -1);
  bool Read(const std::string& input_file, FieldColumnDict& field_columns,
	    AngularCoordinate::Sphere sphere = AngularCoordinate::Equatorial,
	    bool verbose = false, uint8_t theta_column = 0,
	    uint8_t phi_column = 1, int8_t weight_column = -1);

  // Equivalent methods as their namesakes in the BaseMap class.
  virtual void Coverage(PixelVector& superpix,
			uint32_t resolution = HPixResolution);
  bool Covering(Map& stomp_map, uint32_t maximum_pixels);
  virtual double FindUnmaskedFraction(Pixel& pix);
  virtual int8_t FindUnmaskedStatus(Pixel& pix);

  // And if we're not interested in the number of pixels, but want a Map
  // equivalent of the area covered by the nodes in the map.  If the map was
  // built with base level nodes at HPixResolution resolution, the
  // results of Coverage and NodeMap will equivalent, albeit in different
  // functional forms.
  void NodeMap(Map& stomp_map);

  // Some getters and setters for the base level resolution and pixel capacity.
  uint32_t Resolution();
  uint16_t PixelCapacity();
  void SetResolution(uint32_t resolution);
  void SetPixelCapacity(int pixel_capacity);

  // Total number of points in the tree map or total number of points in
  // a given base level node.
  uint32_t NPoints(uint32_t k = MaxPixnum);

  // A variation on the above method, returns the number of points associated
  // with the current map that are also contained in the input pixel.
  uint32_t NPoints(Pixel& pix);

  // If we want to extract a copy of all of the points that have been added
  // to this map, this method allows for that.
  void Points(WAngularVector& w_ang);

  // And an associated method that will extract a copy of the points associated
  // with an input pixel.
  void Points(WAngularVector& w_ang, Pixel& pix);

  // Total weight for all of the points in the tree map or total weight in
  // a given base level node.
  double Weight(uint32_t k = MaxPixnum);

  // Likewise, we can provide a similar method for returning the weight
  // associated with an input pixel.
  double Weight(Pixel& pix);

  // And the equivalent functions for FieldTotals...
  double FieldTotal(const std::string& field_name,
			   uint32_t k = MaxPixnum);
  double FieldTotal(const std::string& field_name, Pixel& pix);
  uint16_t NField();
  bool HasFields();
  void FieldNames(std::vector<std::string>& field_names);

  // Total number of base level nodes.
  uint16_t BaseNodes();

  // Total number of all nodes.
  uint16_t Nodes();

  // We need these methods to comply with the BaseMap signature.
  virtual uint32_t Size();
  virtual double Area();
  void CalculateArea();
  virtual uint32_t MinResolution();
  virtual uint32_t MaxResolution();
  virtual uint8_t MinLevel();
  virtual uint8_t MaxLevel();
  virtual bool Empty();
  virtual void Clear();

 private:
  TreeDict tree_map_;
  FieldDict field_total_;
  uint16_t maximum_points_, nodes_;
  uint32_t point_count_, resolution_;
  double weight_, area_;
  bool modified_;
};

} // end namespace Stomp

#endif
