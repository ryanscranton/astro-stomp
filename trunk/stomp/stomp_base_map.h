// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the abstract BaseMap class that serves as the
// basis for all of the *Map objects.  BaseMap sets out the basic functionality
// that all of the *Map classes need to describe a given region on the sky and
// do some basic internal maintenance.  Additionally, BaseMap provides a
// common set of methods for dividing that area up into nearly equal-area,
// similarly-shaped regions.  This functionality is the basis for calculating
// jack-knife errors for our various statistical analyses.

#ifndef STOMP_BASE_MAP_H
#define STOMP_BASE_MAP_H

#include <stdint.h>
#include <vector>
#include <map>
#include "stomp_core.h"
#include "stomp_pixel.h"

namespace Stomp {

class AngularCoordinate;  // class declaration in stomp_angular_coordinate.h
class Pixel;              // class declaration in stomp_pixel.h
class Section;
class RegionMap;
class BaseMap;

typedef std::map<const uint32_t, int16_t> RegionDict;
typedef RegionDict::iterator RegionIterator;
typedef std::pair<RegionIterator, RegionIterator> RegionPair;

typedef std::map<const int16_t, double> RegionAreaDict;
typedef RegionAreaDict::iterator RegionAreaIterator;
typedef std::pair<RegionAreaIterator, RegionAreaIterator> RegionAreaPair;


class Section {
  // This is barely a class.  Really, it's just a small object that's necessary
  // for constructing the sub-regions in the *Map classes.
 public:
  Section();
  ~Section();
  void SetMinStripe(uint32_t stripe);
  void SetMaxStripe(uint32_t stripe);
  uint32_t MinStripe();
  uint32_t MaxStripe();

 private:
  uint32_t stripe_min_, stripe_max_;
};

class RegionMap {
  // This class provides the functionality for dividing the area subtended by a
  // BaseMap-derived object into roughly equal-area, equal-sized regions.  The
  // class is not intended to be instantiated outside of the BaseMap class.

 public:
  RegionMap();
  virtual ~RegionMap();

  // This method initializes the regions on our current map.  There are two
  // parameters: the resolution to use for breaking up the current map and
  // the number of regions we should divide the current map into.  The higher
  // the resolution value, the more precise our split will be at the expense
  // of more memory.  Resolution values above 2048 are forbidden.  If no
  // resolution value is given, the code will attempt to find a reasonable
  // value based on the requested number of regions.
  //
  // The number of regions can't exceed the number of pixels at the specified
  // resolution for obvious reasons.  Likewise, as the number of regions
  // increases, our ability to make them equal area becomes more constrained.
  // So, don't go crazy here.  The return value is the number of regions that
  // we used in the final splitting, in case the specified number had to be
  // reduced to match the available number of pixels.
  uint16_t InitializeRegions(BaseMap* base_map, uint16_t n_region,
			     uint32_t region_resolution = 0);

  // Alternatively, we could import our region map from another BaseMap.  The
  // return value indicates success or failure.
  //
  // WARNING: Using this method will apply the source_map's regionation
  // to the base_map, regardless of the geometry of the base_map.  This means
  // that, if the two BaseMaps have disjoint areas, then the regionation could
  // have unpredictable results.  Applying a Map's regionation to a ScalarMap
  // or TreeMap derived from that Map or points within the Map will have the
  // expected result.  Other use cases are strictly caveat emptor.
  bool InitializeRegions(BaseMap* base_map, BaseMap& source_map);

  // Once we have the map divided into sub-regions, there are number of things
  // we might do.  The simplest would be to take in an AngularCoordinate object
  // and return the index of the sub-region that contained that point.  If
  // the point is not in any of the regions (and hence, outside of the map),
  // then the return value is -1.
  int16_t FindRegion(AngularCoordinate& ang);

  // And finally, a method for removing the current sub-region setup so that
  // a new version can be imposed on the map.  This method is called before
  // InitializeRegions does anything, so two successive calls to
  // InitializeRegions won't cause problems.
  void ClearRegions();

  // Given a pixel index (the Pixnum method in Pixel), return the corresponding
  // region value.
  int16_t Region(uint32_t region_idx);

  // Given a region index, return a PixelVector corresponding to its area.
  void RegionArea(int16_t region, PixelVector& pix);

  // Given a region index, return the area associated with that region.
  double RegionArea(int16_t region);

  // Some getter methods to describe the state of the RegionMap.
  uint16_t NRegion();
  uint32_t Resolution();
  bool Initialized();

  // Return iterators for the set of RegionMap objects.
  RegionIterator Begin();
  RegionIterator End();

 private:
  RegionDict region_map_;
  RegionAreaDict region_area_;
  uint32_t region_resolution_;
  uint16_t n_region_;
};


class BaseMap {
  // This is the abstract base class that all of the map classes will inherit,
  // prototyping all of the basic map functionality.  In particular, it
  // includes all of the functionality for dividing the map area with the
  // RegionMapper sub-class.  This class should never be instantiated directly.
 public:
  BaseMap();
  virtual ~BaseMap();

  // These four methods are the core methods required for running the
  // RegionMapper code.  It only calls Coverage directly, but the Weight value
  // stored in the Pixels returned indicates the fraction of that pixel
  // contained in the map.  This calculation in turn requires FindUnmaskedArea
  // and FindUnmaskedStatus.  Covering is added as useful description of the
  // overall BaseMap geometry.  See the Map class for more thorough
  // documentation of each of these methods.
  virtual void Coverage(PixelVector& superpix,
			uint32_t resolution = HPixResolution,
			bool calculate_fraction = true);
  virtual double FindUnmaskedFraction(Pixel& pix);
  virtual int8_t FindUnmaskedStatus(Pixel& pix);

  // In addition, each map instance should have some basic book-keeping methods
  // for determining whether the map contains any data (Empty), removing any
  // data (Clear), returning the size of the data load (Size), and giving a
  // notion of the subtended area (Area).
  virtual bool Empty();
  virtual void Clear();
  virtual uint32_t Size();
  virtual double Area();

  // Further, each map needs to give a sense of the range of pixel resolutions
  // involved.  This is important because the RegionMap can't be created at
  // a resolution smaller than the map itself can resolve.
  virtual uint32_t MinResolution();
  virtual uint32_t MaxResolution();
  virtual uint8_t MinLevel();
  virtual uint8_t MaxLevel();

  // These methods all act as wrappers for the RegionMapper object contained
  // in the class.  See that class for documentation.
  uint16_t InitializeRegions(uint16_t n_regions,
			     uint32_t region_resolution = 0);
  bool InitializeRegions(BaseMap& base_map);
  int16_t FindRegion(AngularCoordinate& ang);
  void ClearRegions();
  void RegionArea(int16_t region, PixelVector& pix);
  int16_t Region(uint32_t region_idx);
  double RegionArea(int16_t region);
  uint16_t NRegion();
  uint32_t RegionResolution();
  bool RegionsInitialized();
  RegionIterator RegionBegin();
  RegionIterator RegionEnd();

 private:
  RegionMap region_map_;
};

} // end namespace Stomp

#endif
