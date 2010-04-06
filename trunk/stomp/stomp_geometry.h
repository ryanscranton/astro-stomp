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

namespace Stomp {

class GeometricBound;
class CircleBound;
class WedgeBound;
class PolygonBound;
class Footprint;

typedef std::vector<CircleBound> CircleVector;
typedef CircleVector::iterator CircleIterator;

typedef std::vector<WedgeBound> WedgeVector;
typedef WedgeVector::iterator WedgeIterator;

typedef std::vector<PolygonBound> PolygonVector;
typedef PolygonVector::iterator PolygonIterator;

typedef std::vector<Footprint> FootprintVector;
typedef FootprintVector::iterator FootprintIterator;

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
  // area and bounds.  We need the setters since these values will be stored in
  // the base class.
  void SetArea(double input_area);
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

class CircleBound : public GeometricBound {
  // An example of a derived GeometricBound class.  This implements a simple
  // circular footprint of a given radius (in degrees) around a central
  // angular position.

 public:
  CircleBound(const AngularCoordinate& ang, double radius);
  virtual ~CircleBound();
  virtual bool CheckPoint(AngularCoordinate& ang);
  virtual bool FindAngularBounds();
  virtual bool FindArea();

 private:
  AngularCoordinate ang_;
  double radius_, sin2radius_;
};


class WedgeBound : public GeometricBound {
  // A variant of the CircleBound class.  Instead of pixelizing the entire
  // circle, we only pixelize a wedge from the circle.  The position angle
  // values and radius should be specified in degrees.
 public:
  WedgeBound(const AngularCoordinate& ang, double radius,
	     double position_angle_min, double position_angle_max,
	     AngularCoordinate::Sphere sphere = AngularCoordinate::Survey);
  virtual ~WedgeBound();
  virtual bool CheckPoint(AngularCoordinate& ang);
  virtual bool FindAngularBounds();
  virtual bool FindArea();

 private:
  AngularCoordinate ang_;
  double radius_, sin2radius_, position_angle_min_, position_angle_max_;
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

} // end namespace Stomp

#endif
