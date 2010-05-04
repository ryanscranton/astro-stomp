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

void GeometricBound::SetAngularBounds(double lammin, double lammax,
				      double etamin, double etamax) {
  lammin_ = lammin;
  lammax_ = lammax;
  etamin_ = etamin;
  etamax_ = etamax;
}

void GeometricBound::SetContinuousBounds(bool continuous_bounds) {
  continuous_bounds_ = continuous_bounds;
}

double GeometricBound::Area() {
  return area_;
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

bool GeometricBound::ContinuousBounds() {
  return continuous_bounds_;
}

CircleBound::CircleBound(const AngularCoordinate& ang, double radius) {
  ang_ = ang;
  radius_ = radius;
  sin2radius_ = sin(radius*DegToRad)*sin(radius*DegToRad);

  SetContinuousBounds(true);
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
  bool within_bound = false;

  double costheta =
      ang.UnitSphereX()*ang_.UnitSphereX() +
      ang.UnitSphereY()*ang_.UnitSphereY() +
      ang.UnitSphereZ()*ang_.UnitSphereZ();

  if (DoubleLE(1.0-costheta*costheta, sin2radius_)) within_bound = true;

  return within_bound;
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
  SetContinuousBounds(true);

  FindArea();
  FindAngularBounds();
}

WedgeBound::~WedgeBound() {
  radius_ = sin2radius_ = position_angle_min_ = position_angle_max_ = 0.0;
}

bool WedgeBound::FindAngularBounds() {
  // Start with position directly above center point.
  AngularCoordinate start_ang;

  switch (sphere_) {
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

  double lammin = ang_.Lambda(), lammax = ang_.Lambda();
  double etamin = ang_.Eta(), etamax = ang_.Eta();

  AngularCoordinate rotate_ang;
  start_ang.Rotate(ang_, position_angle_min_, rotate_ang, sphere_);
  if (DoubleLT(rotate_ang.Lambda(), lammin)) lammin = rotate_ang.Lambda();
  if (DoubleGT(rotate_ang.Lambda(), lammax)) lammax = rotate_ang.Lambda();
  if (DoubleLT(rotate_ang.Eta(), etamin)) etamin = rotate_ang.Eta();
  if (DoubleGT(rotate_ang.Eta(), etamax)) etamax = rotate_ang.Eta();

  start_ang.Rotate(ang_, position_angle_max_, rotate_ang, sphere_);
  if (DoubleLT(rotate_ang.Lambda(), lammin)) lammin = rotate_ang.Lambda();
  if (DoubleGT(rotate_ang.Lambda(), lammax)) lammax = rotate_ang.Lambda();
  if (DoubleLT(rotate_ang.Eta(), etamin)) etamin = rotate_ang.Eta();
  if (DoubleGT(rotate_ang.Eta(), etamax)) etamax = rotate_ang.Eta();

  if (DoubleLE(lammin, -90.0)) lammin = -90.0;
  if (DoubleGE(lammax, 90.0)) lammax = 90.0;
  if (DoubleGT(etamin, 180.0)) etamin -= 360.0;
  if (DoubleLT(etamin, -180.0)) etamin += 360.0;

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
  bool within_bound = false;

  double costheta =
    ang.UnitSphereX()*ang_.UnitSphereX() +
    ang.UnitSphereY()*ang_.UnitSphereY() +
    ang.UnitSphereZ()*ang_.UnitSphereZ();

  if (DoubleLE(1.0-costheta*costheta, sin2radius_)) {
    double position_angle = ang_.PositionAngle(ang, sphere_);
    if (DoubleGE(position_angle, position_angle_min_) &&
	DoubleLE(position_angle, position_angle_max_))
      within_bound = true;
  }

  return within_bound;
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

  SetContinuousBounds(true);
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

LongitudeBound::LongitudeBound(double min_longitude, double max_longitude,
			       AngularCoordinate::Sphere sphere) {
  // Cast the input values into AngularCoordinate objects to handle vagaries
  // of signs and longitude bounds.
  AngularCoordinate min_ang, max_ang;
  switch (sphere) {
  case AngularCoordinate::Survey:
    min_ang.SetSurveyCoordinates(0.0, min_longitude);
    max_ang.SetSurveyCoordinates(0.0, max_longitude);

    min_longitude_ = min_ang.Eta();
    max_longitude_ = max_ang.Eta();
    break;
  case AngularCoordinate::Equatorial:
    min_ang.SetEquatorialCoordinates(min_longitude, 0.0);
    max_ang.SetEquatorialCoordinates(max_longitude, 0.0);

    min_longitude_ = min_ang.RA();
    max_longitude_ = max_ang.RA();
    break;
  case AngularCoordinate::Galactic:
    min_ang.SetGalacticCoordinates(min_longitude, 0.0);
    max_ang.SetGalacticCoordinates(max_longitude, 0.0);

    min_longitude_ = min_ang.GalLon();
    max_longitude_ = max_ang.GalLon();
    break;
  }

  if (min_longitude_ < max_longitude_) {
    continuous_longitude_ = true;
  } else {
    continuous_longitude_ = false;
  }

  sphere_ = sphere;

  FindArea();
  FindAngularBounds();
}

LongitudeBound::~LongitudeBound() {
  min_longitude_ = max_longitude_ = 0.0;
}

bool LongitudeBound::CheckPoint(AngularCoordinate& ang) {
  bool inside_bound = false;

  if (continuous_longitude_) {
    switch(sphere_) {
    case AngularCoordinate::Survey:
      if (DoubleLE(ang.Eta(), max_longitude_) &&
	  DoubleGE(ang.Eta(), min_longitude_)) inside_bound = true;
      break;
    case AngularCoordinate::Equatorial:
      if (DoubleLE(ang.RA(), max_longitude_) &&
	  DoubleGE(ang.RA(), min_longitude_)) inside_bound = true;
      break;
    case AngularCoordinate::Galactic:
      if (DoubleLE(ang.GalLon(), max_longitude_) &&
	  DoubleGE(ang.GalLon(), min_longitude_)) inside_bound = true;
      break;
    }
  } else {
    switch(sphere_) {
    case AngularCoordinate::Survey:
      if (DoubleGE(ang.Eta(), min_longitude_) ||
	  DoubleLE(ang.Eta(), max_longitude_)) inside_bound = true;
      break;
    case AngularCoordinate::Equatorial:
      if (DoubleGE(ang.RA(), min_longitude_) ||
	  DoubleLE(ang.RA(), max_longitude_)) inside_bound = true;
      break;
    case AngularCoordinate::Galactic:
      if (DoubleGE(ang.GalLon(), min_longitude_) ||
	  DoubleLE(ang.GalLon(), max_longitude_)) inside_bound = true;
      break;
    }
  }

  return inside_bound;
}

bool LongitudeBound::FindAngularBounds() {
  double lammin = 100.0, lammax = -100.0, etamin = 200.0, etamax = -200.0;
  double lat_min = -90.0;
  double lat_step = 180.0/1000.0;

  switch(sphere_) {
  case AngularCoordinate::Survey:
    lammin = -90.0;
    lammax = 90.0;
    etamin = min_longitude_;
    etamax = max_longitude_;
    break;
  case AngularCoordinate::Equatorial:
    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(min_longitude_, lat_min+lat_step*i,
			    AngularCoordinate::Equatorial);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }

    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(max_longitude_, lat_min+lat_step*i,
			    AngularCoordinate::Equatorial);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }
    break;
  case AngularCoordinate::Galactic:
    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(min_longitude_, lat_min+lat_step*i,
			    AngularCoordinate::Galactic);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }

    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(max_longitude_, lat_min+lat_step*i,
			    AngularCoordinate::Galactic);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }
    break;
  }

  // For large area bounds, the above methods occasionally fail, resulting in
  // bound limits that are too small.  In those cases, we want to just take
  // the entire sphere as our first guess.
  double bound_area = 4.0*Pi*StradToDeg;
  bound_area *= 0.5*(sin(lammin*DegToRad) - sin(lammin*DegToRad));
  bound_area *= (etamax - etamin)/360.0;

  if (bound_area < Area() && Area() > 2.0*Pi) {
    lammin = -90.0;
    lammax = 90.0;
    etamin = -180.0;
    etamax = 180.0;
  }

  SetAngularBounds(lammin,lammax,etamin,etamax);

  return true;
}

bool LongitudeBound::FindArea() {
  double bound_area = 4.0*Pi*StradToDeg;  // start with the whole sphere

  if (continuous_longitude_) {
    bound_area *= (max_longitude_ - max_longitude_)/360.0;
  } else {
    bound_area *= (360.0 - (min_longitude_ - max_longitude_))/360.0;
  }

  SetArea(bound_area);

  return true;
}

double LongitudeBound::LongitudeMin() {
  return min_longitude_;
}

double LongitudeBound::LongitudeMax() {
  return max_longitude_;
}

AngularCoordinate::Sphere LongitudeBound::Sphere() {
  return sphere_;
}

LatitudeBound::LatitudeBound(double min_latitude, double max_latitude,
			       AngularCoordinate::Sphere sphere) {
  min_latitude_ = min_latitude;
  max_latitude_ = max_latitude;

  sphere_ = sphere;

  FindArea();
  FindAngularBounds();
}

LatitudeBound::~LatitudeBound() {
  min_latitude_ = max_latitude_ = 0.0;
}

bool LatitudeBound::CheckPoint(AngularCoordinate& ang) {
  bool inside_bound = false;

  switch(sphere_) {
  case AngularCoordinate::Survey:
    if (DoubleLE(ang.Lambda(), max_latitude_) &&
	DoubleGE(ang.Lambda(), min_latitude_)) inside_bound = true;
    break;
  case AngularCoordinate::Equatorial:
    if (DoubleLE(ang.DEC(), max_latitude_) &&
	DoubleGE(ang.DEC(), min_latitude_)) inside_bound = true;
    break;
  case AngularCoordinate::Galactic:
    if (DoubleLE(ang.GalLat(), max_latitude_) &&
	DoubleGE(ang.GalLat(), min_latitude_)) inside_bound = true;
    break;
  }

  return inside_bound;
}

bool LatitudeBound::FindAngularBounds() {
  double lammin = 100.0, lammax = -100.0, etamin = 200.0, etamax = -200.0;
  double lon_min = 0.0;
  double lon_step = 360.0/1000.0;

  switch(sphere_) {
  case AngularCoordinate::Survey:
    lammin = min_latitude_;
    lammax = max_latitude_;
    etamin = -180.0;
    etamax = 180.0;
    break;
  case AngularCoordinate::Equatorial:
    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(lon_min+lon_step*i, min_latitude_,
			    AngularCoordinate::Equatorial);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }

    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(lon_min+lon_step*i, max_latitude_,
			    AngularCoordinate::Equatorial);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }
    break;
  case AngularCoordinate::Galactic:
    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(lon_min+lon_step*i, min_latitude_,
			    AngularCoordinate::Galactic);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }

    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(lon_min+lon_step*i, max_latitude_,
			    AngularCoordinate::Galactic);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }
    break;
  }

  // For large area bounds, the above methods occasionally fail, resulting in
  // bound limits that are too small.  In those cases, we want to just take
  // the entire sphere as our first guess.
  double bound_area = 4.0*Pi*StradToDeg;
  bound_area *= 0.5*(sin(lammin*DegToRad) - sin(lammin*DegToRad));
  bound_area *= (etamax - etamin)/360.0;

  if (bound_area < Area() && Area() > 2.0*Pi) {
    lammin = -90.0;
    lammax = 90.0;
    etamin = -180.0;
    etamax = 180.0;
  }

  SetAngularBounds(lammin,lammax,etamin,etamax);

  return true;
}

bool LatitudeBound::FindArea() {
  double bound_area = 4.0*Pi*StradToDeg;  // start with the whole sphere

  double z_min = sin(min_latitude_*DegToRad);
  double z_max = sin(max_latitude_*DegToRad);

  bound_area *= (z_max - z_min)/2.0;

  SetArea(bound_area);

  return true;
}

double LatitudeBound::LatitudeMin() {
  return min_latitude_;
}

double LatitudeBound::LatitudeMax() {
  return max_latitude_;
}

AngularCoordinate::Sphere LatitudeBound::Sphere() {
  return sphere_;
}

LatLonBound::LatLonBound(double min_latitude, double max_latitude,
			 double min_longitude, double max_longitude,
			 AngularCoordinate::Sphere sphere) {
  min_latitude_ = min_latitude;
  max_latitude_ = max_latitude;

  // Cast the input values into AngularCoordinate objects to handle vagaries
  // of signs and longitude bounds.
  AngularCoordinate min_ang, max_ang;
  switch (sphere) {
  case AngularCoordinate::Survey:
    min_ang.SetSurveyCoordinates(0.0, min_longitude);
    max_ang.SetSurveyCoordinates(0.0, max_longitude);

    min_longitude_ = min_ang.Eta();
    max_longitude_ = max_ang.Eta();
    break;
  case AngularCoordinate::Equatorial:
    min_ang.SetEquatorialCoordinates(min_longitude, 0.0);
    max_ang.SetEquatorialCoordinates(max_longitude, 0.0);

    min_longitude_ = min_ang.RA();
    max_longitude_ = max_ang.RA();
    break;
  case AngularCoordinate::Galactic:
    min_ang.SetGalacticCoordinates(min_longitude, 0.0);
    max_ang.SetGalacticCoordinates(max_longitude, 0.0);

    min_longitude_ = min_ang.GalLon();
    max_longitude_ = max_ang.GalLon();
    break;
  }

  if (min_longitude_ < max_longitude_) {
    continuous_longitude_ = true;
  } else {
    continuous_longitude_ = false;
  }

  sphere_ = sphere;

  FindArea();
  FindAngularBounds();
}

LatLonBound::~LatLonBound() {
  min_longitude_ = max_longitude_ = 0.0;
  min_latitude_ = max_latitude_ = 0.0;
}

bool LatLonBound::CheckPoint(AngularCoordinate& ang) {
  bool inside_bound = false;

  switch(sphere_) {
  case AngularCoordinate::Survey:
    if (DoubleLE(ang.Lambda(), max_latitude_) &&
	DoubleGE(ang.Lambda(), min_latitude_)) inside_bound = true;
    break;
  case AngularCoordinate::Equatorial:
    if (DoubleLE(ang.DEC(), max_latitude_) &&
	DoubleGE(ang.DEC(), min_latitude_)) inside_bound = true;
    break;
  case AngularCoordinate::Galactic:
    if (DoubleLE(ang.GalLat(), max_latitude_) &&
	DoubleGE(ang.GalLat(), min_latitude_)) inside_bound = true;
    break;
  }

  if (inside_bound) {
    inside_bound = false;

    if (continuous_longitude_) {
      switch(sphere_) {
      case AngularCoordinate::Survey:
	if (DoubleLE(ang.Eta(), max_longitude_) &&
	    DoubleGE(ang.Eta(), min_longitude_)) inside_bound = true;
	break;
      case AngularCoordinate::Equatorial:
	if (DoubleLE(ang.RA(), max_longitude_) &&
	    DoubleGE(ang.RA(), min_longitude_)) inside_bound = true;
	break;
      case AngularCoordinate::Galactic:
	if (DoubleLE(ang.GalLon(), max_longitude_) &&
	    DoubleGE(ang.GalLon(), min_longitude_)) inside_bound = true;
	break;
      }
    } else {
      switch(sphere_) {
      case AngularCoordinate::Survey:
	if (DoubleGE(ang.Eta(), min_longitude_) ||
	    DoubleLE(ang.Eta(), max_longitude_)) inside_bound = true;
	break;
      case AngularCoordinate::Equatorial:
	if (DoubleGE(ang.RA(), min_longitude_) ||
	    DoubleLE(ang.RA(), max_longitude_)) inside_bound = true;
	break;
      case AngularCoordinate::Galactic:
	if (DoubleGE(ang.GalLon(), min_longitude_) ||
	    DoubleLE(ang.GalLon(), max_longitude_)) inside_bound = true;
	break;
      }
    }
  }

  return inside_bound;
}

bool LatLonBound::FindAngularBounds() {
  double lammin = 100.0, lammax = -100.0, etamin = 200.0, etamax = -200.0;
  double lat_min = min_latitude_;
  double lat_step = (max_latitude_ - min_latitude_)/1000.0;
  double lon_min = min_longitude_;
  double lon_step = (max_longitude_ - min_longitude_)/1000.0;
  if (!continuous_longitude_)
    lon_step = (360.0 - (min_longitude_ - max_longitude_))/1000.0;

  switch(sphere_) {
  case AngularCoordinate::Survey:
    lammin = min_latitude_;
    lammax = max_latitude_;
    etamin = min_longitude_;
    etamax = max_longitude_;
    SetContinuousBounds(continuous_longitude_);
    break;
  case AngularCoordinate::Equatorial:
    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(min_longitude_, lat_min+lat_step*i,
			    AngularCoordinate::Equatorial);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }

    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(max_longitude_, lat_min+lat_step*i,
			    AngularCoordinate::Equatorial);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }
    break;

    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(lon_min+lon_step*i, min_latitude_,
			    AngularCoordinate::Equatorial);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }

    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(lon_min+lon_step*i, max_latitude_,
			    AngularCoordinate::Equatorial);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }
    break;
  case AngularCoordinate::Galactic:
    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(min_longitude_, lat_min+lat_step*i,
			    AngularCoordinate::Galactic);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }

    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(max_longitude_, lat_min+lat_step*i,
			    AngularCoordinate::Galactic);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }
    break;

    double lon_min = min_longitude_;
    double lon_step = (max_longitude_ - min_longitude_)/1000.0;
    if (!continuous_longitude_)
      lon_step = (360.0 - (min_longitude_ - max_longitude_))/1000.0;
    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(lon_min+lon_step*i, min_latitude_,
			    AngularCoordinate::Galactic);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }

    for (uint16_t i=0;i<1000;i++) {
      AngularCoordinate ang(lon_min+lon_step*i, max_latitude_,
			    AngularCoordinate::Galactic);
      if (ang.Lambda() < lammin) lammin = ang.Lambda();
      if (ang.Lambda() > lammax) lammax = ang.Lambda();
      if (ang.Eta() < etamin) etamin = ang.Eta();
      if (ang.Eta() > etamax) etamax = ang.Eta();
    }
    break;
  }

  // For large area bounds, the above methods occasionally fail, resulting in
  // bound limits that are too small.  In those cases, we want to just take
  // the entire sphere as our first guess.
  double bound_area = 4.0*Pi*StradToDeg;
  bound_area *= 0.5*(sin(lammin*DegToRad) - sin(lammin*DegToRad));
  bound_area *= (etamax - etamin)/360.0;

  if (bound_area < Area() && Area() > 2.0*Pi) {
    lammin = -90.0;
    lammax = 90.0;
    etamin = -180.0;
    etamax = 180.0;
  }

  SetAngularBounds(lammin,lammax,etamin,etamax);

  return true;
}

bool LatLonBound::FindArea() {
  double bound_area = 4.0*Pi*StradToDeg;  // start with the whole sphere

  double z_min = sin(min_latitude_*DegToRad);
  double z_max = sin(max_latitude_*DegToRad);

  bound_area *= (z_max - z_min)/2.0;

  if (continuous_longitude_) {
    bound_area *= (max_longitude_ - min_longitude_)/360.0;
  } else {
    bound_area *= (360.0 - min_longitude_ + max_longitude_)/360.0;
  }

  SetArea(bound_area);

  return true;
}

double LatLonBound::LongitudeMin() {
  return min_longitude_;
}

double LatLonBound::LongitudeMax() {
  return max_longitude_;
}

double LatLonBound::LatitudeMin() {
  return min_latitude_;
}

double LatLonBound::LatitudeMax() {
  return max_latitude_;
}

AngularCoordinate::Sphere LatLonBound::Sphere() {
  return sphere_;
}

} // end namespace Stomp