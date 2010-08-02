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
#include "stomp_geometry.h"
#include "stomp_pixel.h"

namespace Stomp {

class AngularCoordinate;  // class declaration in stomp_angular_coordinate.h
class Pixel;              // class declaration in stomp_pixel.h
class PixelOrdering;      // class declaration in stomp_pixel.h
class GeometricBound;     // class declaration in stomp_geometry.h
class RegionBound;
class RegionMap;
class BaseMap;

typedef std::vector<RegionBound> RegionBoundVector;
typedef RegionBoundVector::iterator RegionBoundIterator;

typedef std::map<const uint32_t, int16_t> RegionDict;
typedef RegionDict::iterator RegionIterator;
typedef std::pair<RegionIterator, RegionIterator> RegionPair;

typedef std::map<const int16_t, double> RegionAreaDict;
typedef RegionAreaDict::iterator RegionAreaIterator;
typedef std::pair<RegionAreaIterator, RegionAreaIterator> RegionAreaPair;

typedef std::map<const Pixel, bool, PixelOrdering> CoverageDict;
typedef CoverageDict::iterator CoverageIterator;

struct section {
  // Simple structure to hold our section data.  A section is defined by a
  // minimum and maximum Pixel::Stripe value.
  uint32_t min_stripe;
  uint32_t max_stripe;
};

typedef std::vector<section> SectionVector;

class RegionBound {
  // In general, RegionMaps are created in such a way as to make the size and
  // shape of the jack-knife sections as regular as possible.  However, there
  // may be additional geometric information that we want to consider when we
  // split up the area.  A RegionBound allows the user to use a GeometricBound
  // to specify a region of the BaseMap that should be separately considered
  // when dividing up the overall area.  In general, multiple RegionBounds
  // could be specified (with the remaining unclaimed area being parcelled out
  // in the normal way), although users should avoid having overlapping
  // RegionBounds.
 public:
  RegionBound();
  RegionBound(GeometricBound* bound);
  ~RegionBound();

  // If we use the default constructor, we need to be able to set the
  // GeometricBound associated with our RegionBound.
  void SetGeometricBound(GeometricBound* bound);

  // Setter and getter for the number of sub-regions assigned to this
  // RegionBound.
  void SetNRegion(uint16_t n_region);
  uint16_t NRegion();

  // We need to be able to see if a given coverage Pixel is touched by our
  // GeometricBound as well as possibly decide whether a Pixel should go into
  // our RegionBound or another, based on the fraction of the Pixel contained.
  bool CheckPixel(Pixel& pix);
  double ScorePixel(Pixel& pix);

  // Once we know a Pixel covers part of our GeometricBound, we want to include
  // it in our set of Pixels.  However, there's a non-zero possibility that
  // another GeometricBound might contain a larger fraction of the Pixel, so
  // we also need the ability to remove a Pixel from our set.
  bool AddPixel(Pixel& pix);
  bool RemovePixel(Pixel& pix);
  void ClearPixels();

  // Once we have assembled our coverage Pixels, we need a means for extracting
  // a copy of them.  The output PixelVector will be sorted using LocalOrdering.
  void Coverage(PixelVector& pix);
  uint32_t CoveragePixels();

  // Some quick getters to give us access to the area assigned to the
  // RegionBound in terms of the coverage pixels and the GeometricBound area.
  double CoverageArea();
  double BoundArea();


 private:
  GeometricBound* bound_ptr_;
  CoverageDict coverage_pix_;
  double pixel_area_;
  uint16_t n_region_;
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

  // For a finer degree of control, we can also introduce a set of RegionBounds
  // into our calculations.  The idea here is that we may have sub-regions in
  // our BaseMap that we want to treat separately from the rest of the area
  // (the working example here is a survey which has a variable depth over
  // its area).  The RegionBounds are regionated separately from the remainder
  // of the area, ensuring that their boundaries are respected.  With more
  // constraints on the regionation, it is more likely than normal that this
  // way of doing things will result in sub-regions with unequal areas or
  // odd shapes, so be forewarned.
  uint16_t InitializeRegions(BaseMap* base_map,
			     RegionBoundVector& region_bounds,
			     uint16_t n_region,
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

  // Some internal methods we'll use to handle the actual work of regionation.
  //
  // First up, if we aren't given a specific resolution to use for regionation,
  // we need to find one.
  void _FindRegionResolution(BaseMap* base_map, uint16_t n_region,
			     uint32_t region_resolution);

  // Once we have our regionation resolution worked out, the various methods
  // will generate a set of uniform resolution coverage pixels for us to
  // regionate.  The key to making the regions roughly square is to figure out
  // how wide the BaseMap is relative to its length.  To characterize this,
  // we use the Stripe method from the Pixel class.  First we need to find the
  // set of stripes that span our BaseMap.
  void _FindUniqueStripes(PixelVector& coverage_pix,
			  std::vector<uint32_t>& unique_stripes);

  // Given this list of stripes, we can work out where the break points need
  // to be to make the regions roughly square.  This is complicated a bit by
  // the possibility that the BaseMap area isn't simply connected, but this
  // is something we can handle.  This will yield a vector of sections that
  // encode our breakpoints.
  void _FindSections(std::vector<uint32_t>& unique_stripes,
		     double base_map_area, uint16_t n_region,
		     SectionVector& sectionVec);

  // With our sections defined, we can take our vector of coverage
  // pixels and regionate them.
  void _Regionate(PixelVector& coverage_pix, SectionVector& sectionVec,
		  uint16_t n_region, uint16_t starting_region_index = 0);

  // Check that all coverage pixels have a valid region index assigned to
  // them.
  void _VerifyRegionation(uint16_t n_region);

  // Once we have the map divided into sub-regions, there are number of things
  // we might do.  The simplest would be to take in an AngularCoordinate object
  // and return the index of the sub-region that contained that point.  If
  // the point is not in any of the regions (and hence, outside of the map),
  // then the return value is -1.
  int16_t FindRegion(AngularCoordinate& ang);

  // Likewise for an input Pixel (Pixels with resolution lower than our region
  // map will also return -1 even if they are within the BaseMap).
  int16_t FindRegion(Pixel& pix);

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
  uint16_t InitializeRegions(RegionBoundVector& region_bounds,
			     uint16_t n_region,
			     uint32_t region_resolution = 0);
  bool InitializeRegions(BaseMap& base_map);
  int16_t FindRegion(AngularCoordinate& ang);
  int16_t FindRegion(Pixel& pix);
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
