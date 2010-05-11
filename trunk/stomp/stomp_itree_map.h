// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the IndexedTreeMap class.  The core work of
// pair finding and K nearest neighbor searches is done in the
// IndexedTreePixel class.  However, due to the pixelization scheme used in
// STOMP, there is a maximum pixel size that does not span the entire sphere.
// Hence, a vector of IndexedTreePixels is necessary to describe an arbitrary
// collection of points.  IndexedTreeMap manages that vector of
// IndexedTreePixels, adding them as necessary based on the input points.

#ifndef STOMP_ITREE_MAP_H
#define STOMP_ITREE_MAP_H

#include <stdint.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <algorithm>
#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "stomp_itree_pixel.h"
#include "stomp_base_map.h"

namespace Stomp {

class AngularBin;           // class definition in stomp_angular_bin.h
class Map;                  // class definition in stomp_map.h
class IndexedTreePixel;     // class definition in stomp_itree_pixel.h
class IndexedTreeMap;

typedef std::map<const uint32_t, IndexedTreePixel *> ITreeDict;
typedef ITreeDict::iterator ITreeDictIterator;
typedef std::pair<ITreeDictIterator, ITreeDictIterator> ITreeDictPair;

typedef std::vector<IndexedTreeMap> ITreeMapVector;
typedef ITreeMapVector::iterator ITreeMapIterator;
typedef std::pair<ITreeMapIterator, ITreeMapIterator> ITreeMapPair;

class IndexedTreeMap : public BaseMap {
  // A variation on the TreeMap class, but using IndexedTreePixels instead of
  // TreePixels.

 public:
  friend class NearestNeighborPixel;
  // Since the geometry is specified when points are added to the map, the
  // number of parameters for the constructor is potentially very small.  By
  // default the constructor will set things up to have the base resolution of
  // the tree match the resolution of the superpixels and the maximum number of
  // points per node to be 50.  The latter is somewhat arbitrary.
  IndexedTreeMap(uint32_t resolution=HPixResolution,
		 uint16_t maximum_points=50);

  // Alternatively, we can make use of the Read method (see below) to
  // populate our TreeMap during instantion, given an ascii file containing
  // angular coordinate data.
  IndexedTreeMap(
    const std::string& input_file,
    uint32_t resolution=HPixResolution, uint16_t maximum_points=50,
    AngularCoordinate::Sphere sphere = AngularCoordinate::Equatorial,
    bool verbose = false, uint8_t theta_column = 0,
    uint8_t phi_column = 1, int8_t index_column = -1);
  ~IndexedTreeMap();

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

  // In addition to pair finding, we can also use the tree structure we've
  // built to do efficient nearest neighbor searches.  In the general case,
  // we'll be finding the k nearest neighbors of an input point.  The return
  // value is the number of nodes visited during assemblage.
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

  // Analog of the ClosestMatch method in the TreePixel class, where we're only
  // interested in the best match within a given search radius.  The input
  // maximum radius is in degrees and the return boolean indicates whether an
  // acceptable match was found.
  bool ClosestMatch(AngularCoordinate& ang, double max_distance,
		    IndexedAngularCoordinate& match_ang);

  // For the recursion necessary to do the neighbor finding, we use this
  // internal method.
  void _NeighborRecursion(AngularCoordinate& ang,
			  IndexedTreeNeighbor& neighbor);

  // For the recursion necessary to do the closest match finding, we use this
  // internal method.
  void _MatchRecursion(AngularCoordinate& ang,
		       IndexedTreeNeighbor& neighbor);

  // Add a given point on the sphere to the map.
  bool AddPoint(IndexedAngularCoordinate* ang);

  // The default method for adding IndexedAngularCoordinates to the map
  // takes a pointer to the object.  This means that the map now owns that
  // object and it shouldn't be deleted from the heap except by the map.
  // For cases where we want to retain a copy of the point outside of the
  // map, we provide a second method which takes a reference to the object
  // and creates and stores an internal copy.  The input object can thus be
  // modified or deleted without affecting the map.
  bool AddPoint(IndexedAngularCoordinate& i_ang);

  // Complimentary method for specifying a weight separately when adding a
  // point to the pixel.
  bool AddPoint(AngularCoordinate& ang, uint32_t index);

  // Rather than adding points one by one, we can also take an input file and
  // add those points to the tree.  We can do this with and without also adding
  // Field values to each point from the input file.  If the weight column is
  // not specified, then the Weight is set to unity for each point.  As with
  // the ToWAngularVector methods in the IndexedAngularCoordinate class, the
  // returned boolean indicates success or failure in adding all of the points.
  bool Read(const std::string& input_file,
	    AngularCoordinate::Sphere sphere = AngularCoordinate::Equatorial,
	    bool verbose = false, uint8_t theta_column = 0,
	    uint8_t phi_column = 1, int8_t index_column = -1);

  // Equivalent methods as their namesakes in the BaseMap class.
  virtual void Coverage(PixelVector& superpix,
			uint32_t resolution = HPixResolution,
			bool calculate_fraction = true);
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
  void Points(IAngularVector& i_ang);

  // And an associated method that will extract a copy of the points associated
  // with an input pixel.
  void Points(IAngularVector& i_ang, Pixel& pix);

  // Likewise, we can provide a similar method for returning the indices
  // associated with an input pixel.
  void Indices(Pixel& pix, IndexVector& indices);

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
  ITreeDict tree_map_;
  uint16_t maximum_points_, nodes_;
  uint32_t point_count_, resolution_;
  double area_;
  bool modified_;
};

} // end namespace Stomp

#endif
