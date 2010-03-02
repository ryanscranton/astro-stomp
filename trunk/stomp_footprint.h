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

class FootprintBound {
  // This is the base class for generating footprints.  A footprint object is
  // essentially a scaffolding around the Map class that contains the
  // methods necessary for converting some analytic expression of a region on
  // a sphere into the equivalent Map.  All footprints do roughly the same
  // operations to go about this task, but the details differ based on how the
  // analytic decription is implemented.  This is a true abstract class and
  // should never actually be instantiated.  Instead, you should derive classes
  // from this one that replace the virtual methods with ones that are
  // appropriate to your particular footprint geometric description.

 public:
  FootprintBound();
  virtual ~FootprintBound();

  // All footprint derived classes need to replace these virtual methods
  // with ones specific to their respective geometries.  You need a method
  // for saying whether or not a point on the sphere is inside or outside of
  // your area, you need a method to give the footprint object an idea of where
  // to start looking for pixels that might be in your footprint and you need
  // a way to calculate the area of your footprint so that the class can figure
  // out how closely the area of its pixelization of your footprint matches
  // the analytic value.
  virtual bool CheckPoint(AngularCoordinate& ang);
  virtual bool FindAngularBounds();
  virtual bool FindArea();

  // The pixelization method is iteratively adaptive.  First, it tries to find
  // the largest pixels that will likely fit inside the footprint.  Then it
  // checks those pixels against the footprint and keeps the ones that are
  // fully inside the footprint.  For the ones that were at least partially
  // inside the footprint, it refines them to the next resolution level and
  // tests the sub-pixels.  Those that pass are kept, the misses are discarded
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

  // Once we've pixelized the footprint, we want to return a Map
  // representing the results.  This method returns a pointer to that map.
  inline Map::Map* ExportMap() {
    return new Map::Map(pix_);
  }
  inline void ExportMap(Map& stomp_map) {
    stomp_map.Clear();
    stomp_map.Initialize(pix_);
  }
  inline void SetMaxResolution(uint16_t resolution =
			       Stomp::MaxPixelResolution) {
    max_resolution_level_ = Stomp::MostSignificantBit(resolution);
  }

  // Since we store the area and pixelized area in this class, we need methods
  // for setting and getting those values.  Likewise with the weight that will
  // be assigned to the Pixels that will make up the Map that results.
  inline void SetArea(double input_area) {
    area_ = input_area;
  };
  inline double Area() {
    return area_;
  };
  inline void AddToPixelizedArea(uint16_t resolution) {
    pixel_area_ += Pixel::PixelArea(resolution);
  };
  inline double Weight() {
    return weight_;
  };
  inline void SetWeight(double input_weight) {
    weight_ = input_weight;
  };
  inline double PixelizedArea() {
    return pixel_area_;
  };
  inline uint32_t NPixel() {
    return pix_.size();
  };
  inline void SetAngularBounds(double lammin, double lammax,
                               double etamin, double etamax) {
    lammin_ = lammin;
    lammax_ = lammax;
    etamin_ = etamin;
    etamax_ = etamax;
  };
  inline double LambdaMin() {
    return lammin_;
  };
  inline double LambdaMax() {
    return lammax_;
  };
  inline double EtaMin() {
    return etamin_;
  };
  inline double EtaMax() {
    return etamax_;
  };
  inline uint32_t XMin() {
    return x_min_;
  };
  inline uint32_t XMax() {
    return x_max_;
  };
  inline uint32_t YMin() {
    return y_min_;
  };
  inline uint32_t YMax() {
    return y_max_;
  };
  inline PixelIterator Begin() {
    return pix_.begin();
  };
  inline PixelIterator End() {
    return pix_.end();
  };
  inline void Clear() {
    pix_.clear();
  };


 private:
  PixelVector pix_;
  bool found_starting_resolution_, found_xy_bounds_;
  uint8_t max_resolution_level_;
  double area_, pixel_area_, lammin_, lammax_, etamin_, etamax_, weight_;
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
  inline bool DoubleGE(const double x, const double y) {
    static double tolerance = 1.0e-10;
    return (x >= y - tolerance ? true : false);
  };
  inline bool DoubleLE(const double x, const double y) {
    static double tolerance = 1.0e-10;
    return (x <= y + tolerance ? true : false);
  };

 private:
  AngularVector ang_;
  std::vector<double> x_, y_, z_, dot_;
  uint32_t n_vert_;
};

} // end namespace Stomp

#endif
