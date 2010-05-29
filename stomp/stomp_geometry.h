// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// While the general principle of the library is to describe regions on the
// sphere as collections of pixels, it can occasionally be useful to have a
// simple geometric description of a region as well.  The GeometricBound class
// and its derivatives fills this role.  Each GeometricBound instance must be
// able to return its bounded area, its angular bounds (in survey coordinates)
// and indicate whether an input point is inside or outside of its allowed
// area.  For more complicated geometric tasks (like finding the intersection
// of two bounds), the proper procedure would be to create Maps from the
// GeometricBound objects (with the appropriate constructor) and perform those
// operations on their pixelized counterparts.

#ifndef STOMP_GEOMETRY_H
#define STOMP_GEOMETRY_H

#include <stdint.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "stomp_angular_coordinate.h"
#include "stomp_angular_bin.h"
#include "MersenneTwister.h"

namespace Stomp {

class AngularCoordinate;  // class declaration in stomp_angular_coordinate.h
class AngularBin;         // class declaration in stomp_angular_bin.h
class GeometricBound;
class CircleBound;
class AnnulusBound;
class WedgeBound;
class PolygonBound;
class LongitudeBound;
class LatitudeBound;
class LatLonBound;

typedef std::vector<CircleBound> CircleVector;
typedef CircleVector::iterator CircleIterator;

typedef std::vector<AnnulusBound> AnnulusVector;
typedef AnnulusVector::iterator AnnulusIterator;

typedef std::vector<WedgeBound> WedgeVector;
typedef WedgeVector::iterator WedgeIterator;

typedef std::vector<PolygonBound> PolygonVector;
typedef PolygonVector::iterator PolygonIterator;

typedef std::vector<LongitudeBound> LongitudeVector;
typedef LongitudeVector::iterator LongitudeIterator;

typedef std::vector<LatitudeBound> LatitudeVector;
typedef LatitudeVector::iterator LatitudeIterator;

typedef std::vector<LatLonBound> LatLonVector;
typedef LatLonVector::iterator LatLonIterator;

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
  // area and bounds as well as a boolean to indicate whether or not the eta
  // bounds are continuous across the eta discontinuity.  We need the setters
  // since these values will be stored in the base class.
  void SetArea(double input_area);
  void SetAngularBounds(double lammin, double lammax,
			double etamin, double etamax);
  void SetContinuousBounds(bool continuous_bounds);

  double Area();
  double LambdaMin();
  double LambdaMax();
  double EtaMin();
  double EtaMax();
  bool ContinuousBounds();

  // We can add a bit more functionality to the class even with the limited
  // amount of information we have.  In particular, we can generate random
  // points within the bound.  This requires a bit more data, but it will enable
  // us to use monte carlo methods for those cases where that's more efficient
  // than pixelization.
  void GenerateRandomPoint(AngularCoordinate& ang);
  void GenerateRandomPoints(AngularVector& angVec, uint32_t n_rand);

 private:
  MTRand mtrand_;
  double area_, lammin_, lammax_, etamin_, etamax_, z_min_, z_max_;
  bool continuous_bounds_, set_bounds_;
};

class CircleBound : public GeometricBound {
  // An example of a derived GeometricBound class.  This implements a simple
  // circular footprint of a given radius (in degrees) around a central
  // angular position.

 public:
  CircleBound(const AngularCoordinate& center_point, double radius);
  virtual ~CircleBound();
  virtual bool CheckPoint(AngularCoordinate& ang);
  virtual bool FindAngularBounds();
  virtual bool FindArea();

 private:
  AngularCoordinate center_point_;
  double radius_, costhetamin_;
};


class AnnulusBound : public GeometricBound {
  // Similar to the CircleBound, but for an annulus.  In order to keep things
  // simple, we require the circles to share a common center point.

 public:
  AnnulusBound(const AngularCoordinate& center_point, double min_radius,
	       double max_radius);
  AnnulusBound(const AngularCoordinate& center_point,
	       AngularBin& angular_bin);
  virtual ~AnnulusBound();
  virtual bool CheckPoint(AngularCoordinate& ang);
  virtual bool FindAngularBounds();
  virtual bool FindArea();

 private:
  AngularCoordinate center_point_;
  double max_radius_, min_radius_, costhetamin_, costhetamax_;
};


class WedgeBound : public GeometricBound {
  // A variant of the CircleBound class.  Instead of pixelizing the entire
  // circle, we only pixelize a wedge from the circle.  The position angle
  // values and radius should be specified in degrees.
 public:
  WedgeBound(const AngularCoordinate& center_point, double radius,
	     double position_angle_min, double position_angle_max,
	     AngularCoordinate::Sphere sphere = AngularCoordinate::Survey);
  virtual ~WedgeBound();
  virtual bool CheckPoint(AngularCoordinate& ang);
  virtual bool FindAngularBounds();
  virtual bool FindArea();

 private:
  AngularCoordinate center_point_;
  double radius_, costhetamin_, position_angle_min_, position_angle_max_;
  AngularCoordinate::Sphere sphere_;
};


class PolygonBound : public GeometricBound {
  // Another derived GeometricBound class, this one for a spherical polygon
  // represented by a vector of vertices.  In this case, the vertices need to
  // be in clockwise order as seen from outside the sphere, i.e. the order
  // is right-handed when your thumb is pointed towards the center of the
  // sphere.

 public:
  PolygonBound(AngularVector& ang);
  virtual ~PolygonBound();
  virtual bool CheckPoint(AngularCoordinate& ang);
  virtual bool FindAngularBounds();
  virtual bool FindArea();

 private:
  AngularVector ang_;
  std::vector<double> x_, y_, z_, dot_;
  uint32_t n_vert_;
};

class LongitudeBound : public GeometricBound {
  // Simpler than the PolygonBound class, this bound simply includes all of the
  // area on the sphere between a given set of bounds in longitude, for a
  // particular set of spherical coordinates.  Input values of longitude should
  // be in degrees.

 public:
  LongitudeBound(double min_longitude, double max_longitude,
		 AngularCoordinate::Sphere sphere);
  virtual ~LongitudeBound();
  virtual bool CheckPoint(AngularCoordinate& ang);
  virtual bool FindAngularBounds();
  virtual bool FindArea();

  double LongitudeMin();
  double LongitudeMax();
  AngularCoordinate::Sphere Sphere();

 private:
  double min_longitude_, max_longitude_;
  bool continuous_longitude_;
  AngularCoordinate::Sphere sphere_;
};

class LatitudeBound : public GeometricBound {
  // Similar to the LongitudeBound class, but for a given range in latitude.

 public:
  LatitudeBound(double min_latitude, double max_latitude,
		 AngularCoordinate::Sphere sphere);
  virtual ~LatitudeBound();
  virtual bool CheckPoint(AngularCoordinate& ang);
  virtual bool FindAngularBounds();
  virtual bool FindArea();

  double LatitudeMin();
  double LatitudeMax();
  AngularCoordinate::Sphere Sphere();

 private:
  double min_latitude_, max_latitude_;
  AngularCoordinate::Sphere sphere_;
};

class LatLonBound : public GeometricBound {
  // Similar to the LongitudeBound class, but for a given range in both
  // longitude and latitude.

 public:
  LatLonBound(double min_latitude, double max_latitude,
	      double min_longitude, double max_longitude,
	      AngularCoordinate::Sphere sphere);
  virtual ~LatLonBound();
  virtual bool CheckPoint(AngularCoordinate& ang);
  virtual bool FindAngularBounds();
  virtual bool FindArea();

  double LongitudeMin();
  double LongitudeMax();
  double LatitudeMin();
  double LatitudeMax();
  AngularCoordinate::Sphere Sphere();

 private:
  double min_latitude_, max_latitude_, min_longitude_, max_longitude_;
  bool continuous_longitude_;
  AngularCoordinate::Sphere sphere_;
};

} // end namespace Stomp

#endif
