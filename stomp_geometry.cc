// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// The GeometricBound class is intended to be an analytic description of a
// particular region on the sky.  This is in contrast to pixelized *Map classes
// where the area is approximated via some pixelization scheme.  GeometricBound
// is an abstract class intended to represent the minimal functionality required
// for all of the various derived classes in order for them to interface with
// the rest of the library, particularly the Map object.

#include <stdint.h>
#include <iostream>
#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "stomp_geometry.h"

namespace Stomp {

GeometricBound::GeometricBound() {
  FindArea();
  FindAngularBounds();
}

GeometricBound::~GeometricBound() {
  area_ = 0.0;
  lammin_ = etamin_ = 200.0;
  lammax_ = etamax_ = -200.0;
}

bool GeometricBound::CheckPoint(AngularCoordinate& ang) {
  return true;
}

bool GeometricBound::FindAngularBounds() {
  lammin_ = -90.0;
  lammax_ = 90.0;
  etamin_ = -180.0;
  etamax_ = 180.0;

  return true;
}

bool GeometricBound::FindArea() {
  return HPixArea*MaxSuperpixnum;
}

void GeometricBound::SetArea(double input_area) {
  area_ = input_area;
}

double GeometricBound::Area() {
  return area_;
}

void GeometricBound::SetAngularBounds(double lammin, double lammax,
				      double etamin, double etamax) {
  lammin_ = lammin;
  lammax_ = lammax;
  etamin_ = etamin;
  etamax_ = etamax;
}

double GeometricBound::LambdaMin() {
  return lammin_;
}

double GeometricBound::LambdaMax() {
  return lammax_;
}

double GeometricBound::EtaMin() {
  return etamin_;
}

double GeometricBound::EtaMax() {
  return etamax_;
}

CircleBound::CircleBound(const AngularCoordinate& ang, double radius) {
  ang_ = ang;
  radius_ = radius;
  sin2radius_ = sin(radius*DegToRad)*sin(radius*DegToRad);

  FindArea();
  FindAngularBounds();
}

CircleBound::~CircleBound() {
  radius_ = sin2radius_ = 0.0;
}

bool CircleBound::FindAngularBounds() {

  double lammin = ang_.Lambda() - radius_;
  if (DoubleLE(lammin, -90.0)) lammin = -90.0;

  double lammax = ang_.Lambda() + radius_;
  if (DoubleGE(lammax, 90.0)) lammax = 90.0;

  // double eta_multiplier =
  // AngularCoordinate::EtaMultiplier(0.5*(lammax+lammin));
  double eta_multiplier = 1.0;

  double etamin = ang_.Eta() - radius_*eta_multiplier;
  if (DoubleGT(etamin, 180.0)) etamin -= 360.0;
  if (DoubleLT(etamin, -180.0)) etamin += 360.0;

  double etamax = ang_.Eta() + radius_*eta_multiplier;
  if (DoubleGT(etamax, 180.0)) etamax -= 360.0;
  if (DoubleLT(etamax, -180.0)) etamax += 360.0;

  SetAngularBounds(lammin,lammax,etamin,etamax);

  return true;
}

bool CircleBound::FindArea() {
  SetArea((1.0 -
           cos(radius_*DegToRad))*2.0*Pi*StradToDeg);
  return true;
}

bool CircleBound::CheckPoint(AngularCoordinate& ang) {

  double costheta =
      ang.UnitSphereX()*ang_.UnitSphereX() +
      ang.UnitSphereY()*ang_.UnitSphereY() +
      ang.UnitSphereZ()*ang_.UnitSphereZ();

  if (1.0-costheta*costheta <= sin2radius_ + 1.0e-10) return true;

  return false;
}

WedgeBound::WedgeBound(const AngularCoordinate& ang, double radius,
		       double position_angle_min, double position_angle_max,
		       AngularCoordinate::Sphere sphere) {
  ang_ = ang;
  radius_ = radius;
  sin2radius_ = sin(radius*DegToRad)*sin(radius*DegToRad);
  if (DoubleLT(position_angle_min, position_angle_max)) {
    position_angle_max_ = position_angle_max;
    position_angle_min_ = position_angle_min;
  } else {
    position_angle_max_ = position_angle_min;
    position_angle_min_ = position_angle_max;
  }
  sphere_ = sphere;

  FindArea();
  FindAngularBounds();
}

WedgeBound::~WedgeBound() {
  radius_ = sin2radius_ = position_angle_min_ = position_angle_max_ = 0.0;
}

bool WedgeBound::FindAngularBounds() {
  // We should be able to define our bounds based on the three points that
  // define the wedge.
  double lammin = ang_.Lambda();
  double lammax = ang_.Lambda();
  double etamin = ang_.Eta();
  double etamax = ang_.Eta();

  // Now we define a new point directly north of the center of the circle.
  AngularCoordinate start_ang;
  switch(sphere_) {
  case AngularCoordinate::Survey:
    start_ang.SetSurveyCoordinates(ang_.Lambda()+radius_, ang_.Eta());
    break;
  case AngularCoordinate::Equatorial:
    start_ang.SetEquatorialCoordinates(ang_.RA(), ang_.DEC()+radius_);
    break;
  case AngularCoordinate::Galactic:
    start_ang.SetGalacticCoordinates(ang_.GalLon(), ang_.GalLat()+radius_);
    break;
  }

  // Using that point as a reference, we can rotate to our minimum and
  // maximum position angles.
  AngularCoordinate new_ang;
  start_ang.Rotate(ang_, position_angle_min_, new_ang);

  if (DoubleLE(new_ang.Lambda(), lammin)) lammin = new_ang.Lambda();
  if (DoubleGE(new_ang.Lambda(), lammax)) lammax = new_ang.Lambda();

  if (DoubleLE(new_ang.Eta(), etamin)) etamin = new_ang.Eta();
  if (DoubleGE(new_ang.Eta(), etamax)) etamax = new_ang.Eta();

  start_ang.Rotate(ang_, position_angle_max_, new_ang);

  if (DoubleLE(new_ang.Lambda(), lammin)) lammin = new_ang.Lambda();
  if (DoubleGE(new_ang.Lambda(), lammax)) lammax = new_ang.Lambda();

  if (DoubleLE(new_ang.Eta(), etamin)) etamin = new_ang.Eta();
  if (DoubleGE(new_ang.Eta(), etamax)) etamax = new_ang.Eta();

  SetAngularBounds(lammin,lammax,etamin,etamax);

  return true;
}

bool WedgeBound::FindArea() {
  double circle_area =
    (1.0 - cos(radius_*DegToRad))*2.0*Pi*StradToDeg;
  SetArea(circle_area*(position_angle_max_ - position_angle_min_)/360.0);
  return true;
}

bool WedgeBound::CheckPoint(AngularCoordinate& ang) {

  double costheta =
    ang.UnitSphereX()*ang_.UnitSphereX() +
    ang.UnitSphereY()*ang_.UnitSphereY() +
    ang.UnitSphereZ()*ang_.UnitSphereZ();

  if (1.0-costheta*costheta <= sin2radius_ + 1.0e-10) {
    double position_angle = ang_.PositionAngle(ang);
    if (DoubleGE(position_angle, position_angle_min_) &&
	DoubleLE(position_angle, position_angle_max_))
      return true;
  }
  return false;
}

PolygonBound::PolygonBound(AngularVector& ang) {

  for (AngularIterator iter=ang.begin();iter!=ang.end();++iter)
    ang_.push_back(*iter);

  n_vert_ = ang_.size();

  x_.reserve(n_vert_);
  y_.reserve(n_vert_);
  z_.reserve(n_vert_);
  dot_.reserve(n_vert_);

  for (uint32_t i=0;i<n_vert_;i++) {

    std::vector<double> tmp_x, tmp_y, tmp_z;

    for (uint32_t j=0;j<n_vert_;j++) {
      tmp_x.push_back(ang_[j].UnitSphereX());
      tmp_y.push_back(ang_[j].UnitSphereY());
      tmp_z.push_back(ang_[j].UnitSphereZ());
    }

    for (uint32_t j=0;j<n_vert_;j++) {
      if (j == n_vert_ - 1) {
        x_.push_back(tmp_y[j]*tmp_z[0] - tmp_y[0]*tmp_z[j]);
        y_.push_back(tmp_z[j]*tmp_x[0] - tmp_z[0]*tmp_x[j]);
        z_.push_back(tmp_x[j]*tmp_y[0] - tmp_x[0]*tmp_y[j]);
      } else {
        x_.push_back(tmp_y[j]*tmp_z[j+1] - tmp_y[j+1]*tmp_z[j]);
        y_.push_back(tmp_z[j]*tmp_x[j+1] - tmp_z[j+1]*tmp_x[j]);
        z_.push_back(tmp_x[j]*tmp_y[j+1] - tmp_x[j+1]*tmp_y[j]);
      }

      double amplitude = sqrt(x_[j]*x_[j] + y_[j]*y_[j] + z_[j]*z_[j]);

      x_[j] /= amplitude;
      y_[j] /= amplitude;
      z_[j] /= amplitude;

      dot_.push_back(1.0); // This assumes that we're not at constant DEC.
    }
  }

  FindArea();
  FindAngularBounds();
}

PolygonBound::~PolygonBound() {
  ang_.clear();
  x_.clear();
  y_.clear();
  z_.clear();
  dot_.clear();
  n_vert_ = 0;
}

bool PolygonBound::FindAngularBounds() {

  double lammin = 100.0, lammax = -100.0, etamin = 200.0, etamax = -200.0;

  for (uint32_t i=0;i<n_vert_;i++) {
    if (ang_[i].Lambda() < lammin) lammin = ang_[i].Lambda();
    if (ang_[i].Lambda() > lammax) lammax = ang_[i].Lambda();
    if (ang_[i].Eta() < etamin) etamin = ang_[i].Eta();
    if (ang_[i].Eta() > etamax) etamax = ang_[i].Eta();
  }

  SetAngularBounds(lammin,lammax,etamin,etamax);

  return true;
}

bool PolygonBound::FindArea() {

  double sum = 0.0;

  for (uint32_t j=0,k=1;j<n_vert_;j++,k++) {
    if (k == n_vert_) k = 0;

    double cm = (-x_[j]*x_[k] - y_[j]*y_[k] - z_[j]*z_[k]);

    sum += acos(cm);
  }

  double tmp_area = (sum - (n_vert_ - 2)*Pi)*StradToDeg;

  if (tmp_area > 4.0*Pi*StradToDeg) {
    std::cout << "Polygon area is over half the sphere.  This is bad.\n";
    return false;
  }

  SetArea(tmp_area);

  return true;
}

bool PolygonBound::CheckPoint(AngularCoordinate& ang) {
  bool in_polygon = true;

  uint32_t n=0;
  while ((n < n_vert_) && (in_polygon)) {

    in_polygon = false;
    double dot = 1.0 - x_[n]*ang.UnitSphereX() -
        y_[n]*ang.UnitSphereY() - z_[n]*ang.UnitSphereZ();
    if (DoubleLE(dot_[n],0.0)) {
      if (DoubleLE(fabs(dot_[n]), dot)) in_polygon = true;
    } else {
      if (DoubleGE(dot_[n], dot)) in_polygon = true;
    }
    n++;
  }

  return in_polygon;
}

} // end namespace Stomp
