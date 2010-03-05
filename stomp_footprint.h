// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
//  Stuff about footprints

#ifndef STOMP_FOOTPRINT_H
#define STOMP_FOOTPRINT_H

#include <stdint.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "stomp_angular_coordinate.h"
#include "stomp_pixel.h"

namespace Stomp {

class Map;              // class definition in stomp_map.h
class FootprintBound;
class CircleBound;
class WedgeBound;
class PolygonBound;

typedef std::vector<CircleBound> CircleVector;
typedef CircleVector::iterator CircleIterator;

typedef std::vector<WedgeBound> WedgeVector;
typedef WedgeVector::iterator WedgeIterator;

typedef std::vector<PolygonBound> PolygonVector;
typedef PolygonVector::iterator PolygonIterator;

class GeometricBound {
  // The starting point for a footprint is some analytic geometric bound.  This
  // object needs to be able to do three things:
  //
  // * Given a point on the sky, return a boolean indicating whether that
  //   point is inside or outside the bound.
  // * Determine its angular extent in Survey coordinates.
  // * Determine its area.
  //
  // In principle, GeometricBounds could be complex objects composed of
  // multiple GeometricBounds which could be combined to include or exclude
  // regions.  However, the difficulty in keeping track of the bounds and area
  // of such an object will limit our ambitions for the moment.
  //
  // Either way, the GeometricBound class is intended to be an abstract class.
  // The derived classes will handle the actual implementation of these methods.
 public:
  GeometricBound();
  virtual ~GeometricBound();

  // The three methods that all derived classes need to implement.
  virtual bool CheckPoint(AngularCoordinate& ang);
  virtual bool FindAngularBounds();
  virtual bool FindArea();

  // And some simple getters and setters for the values that determine the
  // area and bounds.
  void SetArea(); // need this since we store in the area in the base class.
  double Area();
  void SetAngularBounds(double lammin, double lammax,
			double etamin, double etamax);
  double LambdaMin();
  double LambdaMax();
  double EtaMin();
  double EtaMax();

 private:
  double area_, lammin_, lammax_, etamin_, etamax_;
};


class Footprint {
  // This is the base class for generating footprints.  A footprint object is
  // essentially a scaffolding around the Map class that contains the
  // methods necessary for converting some analytic expression of a region on
  // a sphere into the equivalent Map.  Once this is done this Map can be
  // exported from the Footprint object or queried from within.
 public:
  FootprintBound();

  // In addition to the default constructor, we can also instantiate it with
  // the GeometricBound object we want pixelized and initiate pixelizaiton
  // right away with the input weight value.
  FootprintBound(GeometricBound& geometric_bound, double weight = 1.0,
		 bool pixelize_bound = true);
  virtual ~FootprintBound();

  // If we didn't instantiate with one, we need to set the GeometricBound that
  // we're going to pixelize.
  void SetBound(GeometricBound& geometric_bound);

  // Setter and getter for the weight assigned to the pixels in the Map.
  void SetWeight(double weight);
  double Weight();

  // The pixelization method is iteratively adaptive.  First, it tries to find
  // the largest pixels that will likely fit inside the footprint.  Then it
  // checks those pixels against the footprint and keeps the ones that are
  // fully inside the footprint.  For the ones that were at least partially
  // inside the footprint, it refines them to the next resolution level and
  // tests their sub-pixels.  Those that pass are kept, the misses are discarded
  // and the partials are refined again.  This continues until we reach the
  // maximum resolution level, at which point we keep enough of the partials
  // to match the footprint's area.  When doing the job of pixelizing a given
  // footprint, these three methods should be called subsquently, with the
  // output of FindStartingResolution fed into FindXYBounds as its
  // argument.  A false return value for either FindXYBounds or Pixelize
  // indicates a failure in the corresponding step.
  //
  // Alternatively, just call Pixelize() and it will call the other two
  // routines as necessary.
  uint8_t FindStartingResolutionLevel();
  bool FindXYBounds(const uint8_t resolution_level);
  bool Pixelize();

  // Part of the pixelization process is figuring out what fraction of a
  // given pixel is within the bounds delineated by the footprint's geometry.
  // Pixels are scored on a scale from -1 <= score <= 0, with -1 indicating
  // that the pixel is completely inside of the bounds and 0 indicating that
  // it's completely outside.  This allows one to sort pixels by their score
  // and keep the ones that are closest to being completely inside the
  // footprint bounds.
  double ScorePixel(Pixel& pix);

  // Since we're translating an analytic bound into a set of pixels, we need
  // to choose the level of fidelity we want to use.  A higher resolution means
  // a better approximation of the analytic bound, but it also means a larger
  // set of pixels and a longer time to pixelize.
  void SetMaxResolution(uint16_t resolution = MaxPixelResolution);

  // Once we've pixelized the footprint, we want to return a Map
  // representing the results.  This method returns a pointer to that map.
  Map::Map* ExportMap();
  void ExportMap(Map& stomp_map);

  // Alternatively, we can iterate through the Map stored internally as we
  // do with the Map object.
  MapIterator Begin();
  MapIterator End();
  void Iterate(MapIterator* iter);

  // By default we return the area given by the geometric bound.  Alternatively,
  // we can return the area of the pixelized Map.
  double Area();
  double PixelizedArea();

  // Some of the standard methods from the Map class.  We omit any of the
  // methods that would change the geometry of the underlying Map since we want
  // it to match the input GeometricBound.
  bool FindLocation(AngularCoordinate& ang);
  bool FindLocation(AngularCoordinate& ang, double& weight);
  double FindLocationWeight(AngularCoordinate& ang);

  bool Contains(Pixel& pix);
  bool Contains(Map& stomp_map);

  double FindUnmaskedFraction(Pixel& pix);
  void FindUnmaskedFraction(PixelVector& pix,
                            std::vector<double>& unmasked_fraction);
  void FindUnmaskedFraction(PixelVector& pix);
  double FindUnmaskedFraction(Map& stomp_map);

  int8_t FindUnmaskedStatus(Pixel& pix);
  void FindUnmaskedStatus(PixelVector& pix,
                          std::vector<int8_t>& unmasked_status);
  int8_t FindUnmaskedStatus(Map& stomp_map);

  double FindAverageWeight(Pixel& pix);
  void FindAverageWeight(PixelVector& pix,
                         std::vector<double>& average_weight);
  void FindAverageWeight(PixelVector& pix);
  double AverageWeight();

  void FindMatchingPixels(Pixel& pix,
                          PixelVector& match_pix,
                          bool use_local_weights = false);
  void FindMatchingPixels(PixelVector& pix,
                          PixelVector& match_pix,
                          bool use_local_weights = false);

  void Coverage(PixelVector& superpix,
		uint16_t resolution = Stomp::HPixResolution);
  bool Covering(Map& stomp_map, uint32_t maximum_pixels);

  uint16_t InitializeRegions(uint16_t n_regions,
                             uint16_t region_resolution = 0);
  bool InitializeRegions(BaseMap& base_map);
  int16_t FindRegion(AngularCoordinate& ang);
  void ClearRegions();
  void RegionArea(int16_t region, PixelVector& pix);
  int16_t Region(uint32_t region_idx);
  double RegionArea(int16_t region);
  uint16_t NRegion();
  uint16_t RegionResolution();
  bool RegionsInitialized();
  RegionIterator RegionBegin();
  RegionIterator RegionEnd();
  bool RegionOnlyMap(int16_t region_index, Map& stomp_map);
  bool RegionExcludedMap(int16_t region_index, Map& stomp_map);

  void GenerateRandomPoints(AngularVector& ang, uint32_t n_point = 1,
                            bool use_weighted_sampling = false);
  void GenerateRandomPoints(WAngularVector& ang, WAngularVector& input_ang);
  void GenerateRandomPoints(WAngularVector& ang, std::vector<double>& weights);

  bool Write(std::string& OutputFile, bool hpixel_format = true,
             bool weighted_map = true);

  void ScaleWeight(const double weight_scale);
  void AddConstantWeight(const double add_weight);
  void InvertWeight();

  void Pixels(PixelVector& pix, uint32_t superpixnum = Stomp::MaxSuperpixnum);

  bool ContainsSuperpixel(uint32_t superpixnum);

  uint16_t MinResolution();
  uint16_t MinResolution(uint32_t superpixnum);
  uint16_t MaxResolution();
  uint16_t MaxResolution(uint32_t superpixnum);
  double MinWeight();
  double MinWeight(uint32_t superpixnum);
  double MaxWeight();
  double MaxWeight(uint32_t superpixnum);
  uint32_t Size();
  uint32_t Size(uint32_t superpixnum);
  bool Empty();
  uint32_t PixelCount(uint16_t resolution);

  // Reset everything and delete internal Map and GeometricBound.
  void Clear();

 private:
  Map stomp_map_;
  GeometricBound bound_;
  bool found_starting_resolution_, found_xy_bounds_;
  uint8_t max_resolution_level_;
  double weight_;
  uint32_t x_min_, x_max_, y_min_, y_max_;
};


class CircleBound : public FootprintBound {
  // An example of a derived FootprintBound class.  This implements a simple
  // circular footprint of a given radius (in degrees) around a central
  // angular position.

 public:
  CircleBound(const AngularCoordinate& ang, double radius, double weight);
  virtual ~CircleBound();
  virtual bool CheckPoint(AngularCoordinate& ang);
  virtual bool FindAngularBounds();
  virtual bool FindArea();

 private:
  AngularCoordinate ang_;
  double radius_, sin2radius_;
};


class WedgeBound : public FootprintBound {
  // A variant of the CircleBound class.  Instead of pixelizing the entire
  // circle, we only pixelize a wedge from the circle.  The position angle
  // values and radius should be specified in degrees.
 public:
  WedgeBound(const AngularCoordinate& ang, double radius,
	     double position_angle_min, double position_angle_max,
	     double weight, AngularCoordinate::Sphere sphere =
	     AngularCoordinate::Survey);
  virtual ~WedgeBound();
  virtual bool CheckPoint(AngularCoordinate& ang);
  virtual bool FindAngularBounds();
  virtual bool FindArea();

 private:
  AngularCoordinate ang_;
  double radius_, sin2radius_, position_angle_min_, position_angle_max_;
  AngularCoordinate::Sphere sphere_;
};


class PolygonBound : public FootprintBound {
  // Another derived FootprintBoundClass, this one for a spherical polygon
  // represented by a vector of vertices.  In this case, the vertices need to
  // be in clockwise order as seen from outside the sphere, i.e. the order
  // is right-handed when your thumb is pointed towards the center of the
  // sphere.

 public:
  PolygonBound(AngularVector& ang, double weight);
  virtual ~PolygonBound();
  virtual bool CheckPoint(AngularCoordinate& ang);
  virtual bool FindAngularBounds();
  virtual bool FindArea();

 private:
  AngularVector ang_;
  std::vector<double> x_, y_, z_, dot_;
  uint32_t n_vert_;
};

} // end namespace Stomp

#endif
