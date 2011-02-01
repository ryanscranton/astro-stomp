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
#include "stomp_pixel.h"
#include "stomp_angular_coordinate.h"
#include "stomp_geometry.h"

namespace Stomp {

GeometricBound::GeometricBound() {
  set_bounds_ = false;
  mtrand_.seed();
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
  SetAngularBounds(-90.0, 90.0, -180.0, 180.0);
  return true;
}

bool GeometricBound::FindArea() {
  return HPixArea*MaxSuperpixnum;
}

bool GeometricBound::CheckPixel(Pixel& pix) {
  double inv_nx = 1.0/static_cast<double>(Nx0*pix.Resolution());
  double inv_ny = 1.0/static_cast<double>(Ny0*pix.Resolution());
  double x = static_cast<double>(pix.PixelX());
  double y = static_cast<double>(pix.PixelY());

  double lammid = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.5)*inv_ny);
  double lammin = 90.0 - RadToDeg*acos(1.0-2.0*(y+1.0)*inv_ny);
  double lammax = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.0)*inv_ny);
  double lam_quart = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.75)*inv_ny);
  double lam_three = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.25)*inv_ny);

  double etamid = RadToDeg*(2.0*Pi*(x+0.5))*inv_nx + EtaOffSet;
  if (DoubleGE(etamid, 180.0)) etamid -= 360.0;
  if (DoubleLE(etamid, -180.0)) etamid += 360.0;

  double etamin = RadToDeg*(2.0*Pi*(x+0.0))*inv_nx + EtaOffSet;
  if (DoubleGE(etamin, 180.0)) etamin -= 360.0;
  if (DoubleLE(etamin, -180.0)) etamin += 360.0;

  double etamax = RadToDeg*(2.0*Pi*(x+1.0))*inv_nx + EtaOffSet;
  if (DoubleGE(etamax, 180.0)) etamax -= 360.0;
  if (DoubleLE(etamax, -180.0)) etamax += 360.0;

  double eta_quart = RadToDeg*(2.0*Pi*(x+0.25))*inv_nx + EtaOffSet;
  if (DoubleGE(eta_quart, 180.0)) eta_quart -= 360.0;
  if (DoubleLE(eta_quart, -180.0)) eta_quart += 360.0;

  double eta_three = RadToDeg*(2.0*Pi*(x+0.75))*inv_nx + EtaOffSet;
  if (DoubleGE(eta_three, 180.0)) eta_three -= 360.0;
  if (DoubleLE(eta_three, -180.0)) eta_three += 360.0;

  bool within_bound = false;

  AngularCoordinate ang(lammid, etamid, AngularCoordinate::Survey);
  if (CheckPoint(ang)) within_bound = true;

  if (!within_bound) {
    ang.SetSurveyCoordinates(lam_quart,etamid);
    if (CheckPoint(ang)) within_bound = true;
    ang.SetSurveyCoordinates(lam_three,etamid);
    if (CheckPoint(ang)) within_bound = true;
    ang.SetSurveyCoordinates(lammid,eta_quart);
    if (CheckPoint(ang)) within_bound = true;
    ang.SetSurveyCoordinates(lammid,eta_quart);
    if (CheckPoint(ang)) within_bound = true;
  }

  if (!within_bound) {
    ang.SetSurveyCoordinates(lam_quart,eta_quart);
    if (CheckPoint(ang)) within_bound = true;
    ang.SetSurveyCoordinates(lam_three,eta_quart);
    if (CheckPoint(ang)) within_bound = true;
    ang.SetSurveyCoordinates(lam_quart,eta_three);
    if (CheckPoint(ang)) within_bound = true;
    ang.SetSurveyCoordinates(lam_three,eta_three);
    if (CheckPoint(ang)) within_bound = true;
  }

  if (!within_bound) {
    ang.SetSurveyCoordinates(lammid,etamax);
    if (CheckPoint(ang)) within_bound = true;
    ang.SetSurveyCoordinates(lammid,etamin);
    if (CheckPoint(ang)) within_bound = true;
    ang.SetSurveyCoordinates(lammax,etamid);
    if (CheckPoint(ang)) within_bound = true;
    ang.SetSurveyCoordinates(lammin,etamid);
    if (CheckPoint(ang)) within_bound = true;
  }

  if (!within_bound) {
    ang.SetSurveyCoordinates(lammax,etamax);
    if (CheckPoint(ang)) within_bound = true;
    ang.SetSurveyCoordinates(lammax,etamin);
    if (CheckPoint(ang)) within_bound = true;
    ang.SetSurveyCoordinates(lammin,etamax);
    if (CheckPoint(ang)) within_bound = true;
    ang.SetSurveyCoordinates(lammin,etamin);
    if (CheckPoint(ang)) within_bound = true;
  }

  return within_bound;
}

double GeometricBound::ScorePixel(Pixel& pix) {
  double inv_nx = 1.0/static_cast<double>(Nx0*pix.Resolution());
  double inv_ny = 1.0/static_cast<double>(Ny0*pix.Resolution());
  double x = static_cast<double>(pix.PixelX());
  double y = static_cast<double>(pix.PixelY());

  double lammid = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.5)*inv_ny);
  double lammin = 90.0 - RadToDeg*acos(1.0-2.0*(y+1.0)*inv_ny);
  double lammax = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.0)*inv_ny);
  double lam_quart = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.75)*inv_ny);
  double lam_three = 90.0 - RadToDeg*acos(1.0-2.0*(y+0.25)*inv_ny);

  double etamid = RadToDeg*(2.0*Pi*(x+0.5))*inv_nx + EtaOffSet;
  if (DoubleGE(etamid, 180.0)) etamid -= 360.0;
  if (DoubleLE(etamid, -180.0)) etamid += 360.0;

  double etamin = RadToDeg*(2.0*Pi*(x+0.0))*inv_nx + EtaOffSet;
  if (DoubleGE(etamin, 180.0)) etamin -= 360.0;
  if (DoubleLE(etamin, -180.0)) etamin += 360.0;

  double etamax = RadToDeg*(2.0*Pi*(x+1.0))*inv_nx + EtaOffSet;
  if (DoubleGE(etamax, 180.0)) etamax -= 360.0;
  if (DoubleLE(etamax, -180.0)) etamax += 360.0;

  double eta_quart = RadToDeg*(2.0*Pi*(x+0.25))*inv_nx + EtaOffSet;
  if (DoubleGE(eta_quart, 180.0)) eta_quart -= 360.0;
  if (DoubleLE(eta_quart, -180.0)) eta_quart += 360.0;

  double eta_three = RadToDeg*(2.0*Pi*(x+0.75))*inv_nx + EtaOffSet;
  if (DoubleGE(eta_three, 180.0)) eta_three -= 360.0;
  if (DoubleLE(eta_three, -180.0)) eta_three += 360.0;

  double score = 0.0;

  AngularCoordinate ang(lammid, etamid, AngularCoordinate::Survey);
  if (CheckPoint(ang)) score -= 4.0;

  ang.SetSurveyCoordinates(lam_quart,etamid);
  if (CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lam_three,etamid);
  if (CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lammid,eta_quart);
  if (CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lammid,eta_quart);
  if (CheckPoint(ang)) score -= 3.0;

  ang.SetSurveyCoordinates(lam_quart,eta_quart);
  if (CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lam_three,eta_quart);
  if (CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lam_quart,eta_three);
  if (CheckPoint(ang)) score -= 3.0;
  ang.SetSurveyCoordinates(lam_three,eta_three);
  if (CheckPoint(ang)) score -= 3.0;

  ang.SetSurveyCoordinates(lammid,etamax);
  if (CheckPoint(ang)) score -= 2.0;
  ang.SetSurveyCoordinates(lammid,etamin);
  if (CheckPoint(ang)) score -= 2.0;
  ang.SetSurveyCoordinates(lammax,etamid);
  if (CheckPoint(ang)) score -= 2.0;
  ang.SetSurveyCoordinates(lammin,etamid);
  if (CheckPoint(ang)) score -= 2.0;

  ang.SetSurveyCoordinates(lammax,etamax);
  if (CheckPoint(ang)) score -= 1.0;
  ang.SetSurveyCoordinates(lammax,etamin);
  if (CheckPoint(ang)) score -= 1.0;
  ang.SetSurveyCoordinates(lammin,etamax);
  if (CheckPoint(ang)) score -= 1.0;
  ang.SetSurveyCoordinates(lammin,etamin);
  if (CheckPoint(ang)) score -= 1.0;

  return score/40.0;
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

  z_min_ = sin(DegToRad*lammin_);
  z_max_ = sin(DegToRad*lammax_);
  set_bounds_ = true;
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

void GeometricBound::GenerateRandomPoint(AngularCoordinate& ang) {
  bool keep = false;

  while (!keep) {
    double z = z_min_ + mtrand_.rand(z_max_ - z_min_);
    double lambda = asin(z)*RadToDeg;
    double eta = etamin_ + mtrand_.rand(etamax_ - etamin_);
    ang.SetSurveyCoordinates(lambda, eta);

    if (CheckPoint(ang)) keep = true;
  }

}

void GeometricBound::GenerateRandomPoints(AngularVector& angVec,
					  uint32_t n_rand) {
  if (!angVec.empty()) angVec.clear();
  angVec.reserve(n_rand);

  bool keep = false;
  AngularCoordinate tmp_ang(0.0,0.0);

  for (uint32_t i=0;i<n_rand;i++) {
    while (!keep) {
      double z = z_min_ + mtrand_.rand(z_max_ - z_min_);
      double lambda = asin(z)*RadToDeg;
      double eta = etamin_ + mtrand_.rand(etamax_ - etamin_);

      tmp_ang.SetSurveyCoordinates(lambda, eta);
      if (CheckPoint(tmp_ang)) keep = true;
    }
    angVec.push_back(tmp_ang);
  }
}

CircleBound::CircleBound(const AngularCoordinate& center_point,
			 double radius) {
  center_point_ = center_point;
  radius_ = radius;
  costhetamin_ = cos(radius*DegToRad);

  SetContinuousBounds(true);
  FindArea();
  FindAngularBounds();
}

CircleBound::~CircleBound() {
  radius_ = costhetamin_ = 0.0;
}

bool CircleBound::FindAngularBounds() {

  double lammin = center_point_.Lambda() - radius_;
  if (DoubleLE(lammin, -90.0)) lammin = -90.0;

  double lammax = center_point_.Lambda() + radius_;
  if (DoubleGE(lammax, 90.0)) lammax = 90.0;

  // double eta_multiplier =
  // AngularCoordinate::EtaMultiplier(0.5*(lammax+lammin));
  double eta_multiplier = 1.0;

  double etamin = center_point_.Eta() - radius_*eta_multiplier;
  if (DoubleGT(etamin, 180.0)) etamin -= 360.0;
  if (DoubleLT(etamin, -180.0)) etamin += 360.0;

  double etamax = center_point_.Eta() + radius_*eta_multiplier;
  if (DoubleGT(etamax, 180.0)) etamax -= 360.0;
  if (DoubleLT(etamax, -180.0)) etamax += 360.0;

  SetAngularBounds(lammin,lammax,etamin,etamax);

  return true;
}

bool CircleBound::FindArea() {
  SetArea((1.0 - costhetamin_)*2.0*Pi*StradToDeg);
  return true;
}

bool CircleBound::CheckPoint(AngularCoordinate& ang) {
  return (DoubleGE(center_point_.DotProduct(ang), costhetamin_) ? true : false);
}

AnnulusBound::AnnulusBound(const AngularCoordinate& center_point,
			   double min_radius, double max_radius) {
  center_point_ = center_point;

  if (min_radius > max_radius) {
    min_radius_ = max_radius;
    max_radius_ = min_radius;
  } else {
    min_radius_ = min_radius;
    max_radius_ = max_radius;
  }

  costhetamax_ = cos(min_radius_*DegToRad);
  costhetamin_ = cos(max_radius_*DegToRad);

  SetContinuousBounds(true);
  FindArea();
  FindAngularBounds();
}

AnnulusBound::AnnulusBound(const AngularCoordinate& center_point,
			   AngularBin& angular_bin) {
  center_point_ = center_point;

  min_radius_ = angular_bin.ThetaMin();
  max_radius_ = angular_bin.ThetaMax();

  costhetamax_ = cos(min_radius_*DegToRad);
  costhetamin_ = cos(max_radius_*DegToRad);

  SetContinuousBounds(true);
  FindArea();
  FindAngularBounds();
}

AnnulusBound::~AnnulusBound() {
  min_radius_ = max_radius_ = 0.0;
  costhetamin_ = costhetamax_ = 1.0;
}

bool AnnulusBound::FindAngularBounds() {

  double lammin = center_point_.Lambda() - max_radius_;
  if (DoubleLE(lammin, -90.0)) lammin = -90.0;

  double lammax = center_point_.Lambda() + max_radius_;
  if (DoubleGE(lammax, 90.0)) lammax = 90.0;

  // double eta_multiplier =
  // AngularCoordinate::EtaMultiplier(0.5*(lammax+lammin));
  double eta_multiplier = 1.0;

  double etamin = center_point_.Eta() - max_radius_*eta_multiplier;
  if (DoubleGT(etamin, 180.0)) etamin -= 360.0;
  if (DoubleLT(etamin, -180.0)) etamin += 360.0;

  double etamax = center_point_.Eta() + max_radius_*eta_multiplier;
  if (DoubleGT(etamax, 180.0)) etamax -= 360.0;
  if (DoubleLT(etamax, -180.0)) etamax += 360.0;

  SetAngularBounds(lammin,lammax,etamin,etamax);

  return true;
}

bool AnnulusBound::FindArea() {
  SetArea((costhetamax_ - costhetamin_)*2.0*Pi*StradToDeg);
  return true;
}

bool AnnulusBound::CheckPoint(AngularCoordinate& ang) {
  return (DoubleLE(center_point_.DotProduct(ang), costhetamax_) &&
	  DoubleGE(center_point_.DotProduct(ang), costhetamin_) ? true : false);
}

WedgeBound::WedgeBound(const AngularCoordinate& center_point, double radius,
		       double position_angle_min, double position_angle_max,
		       AngularCoordinate::Sphere sphere) {
  center_point_ = center_point;
  radius_ = radius;
  costhetamin_ = cos(radius*DegToRad);
  if (DoubleLT(position_angle_min, position_angle_max)) {
    position_angle_max_ = position_angle_max;
    position_angle_min_ = position_angle_min;
  } else {
    position_angle_max_ = position_angle_min;
    position_angle_min_ = position_angle_max;
  }
  sphere_ = sphere;

  AngularCoordinate start_ang;
  switch (sphere_) {
  case AngularCoordinate::Survey:
    start_ang.SetSurveyCoordinates(center_point_.Lambda()+radius_,
				   center_point_.Eta());
    break;
  case AngularCoordinate::Equatorial:
    start_ang.SetEquatorialCoordinates(center_point_.RA(),
				       center_point_.DEC()+radius_);
    break;
  case AngularCoordinate::Galactic:
    start_ang.SetGalacticCoordinates(center_point_.GalLon(),
				     center_point_.GalLat()+radius_);
    break;
  }

  AngularCoordinate rotate_ang;
  start_ang.Rotate(center_point_, position_angle_min_, rotate_ang, sphere_);
  double cosphi = center_point_.CosPositionAngle(rotate_ang, sphere_);
  double sinphi = center_point_.SinPositionAngle(rotate_ang, sphere_);
  cosphi_min_ = cosphi_max_ = cosphi;
  sinphi_min_ = sinphi_max_ = sinphi;

  start_ang.Rotate(center_point_, position_angle_max_, rotate_ang, sphere_);
  cosphi = center_point_.CosPositionAngle(rotate_ang, sphere_);
  sinphi = center_point_.SinPositionAngle(rotate_ang, sphere_);
  if (DoubleLT(cosphi, cosphi_min_)) cosphi_min_ = cosphi;
  if (DoubleGT(cosphi, cosphi_max_)) cosphi_max_ = cosphi;
  if (DoubleLT(sinphi, sinphi_min_)) sinphi_min_ = sinphi;
  if (DoubleGT(sinphi, sinphi_max_)) sinphi_max_ = sinphi;

  SetContinuousBounds(true);

  FindArea();
  FindAngularBounds();
}

WedgeBound::~WedgeBound() {
  radius_ = position_angle_min_ = position_angle_max_ = 0.0;
  costhetamin_ = 1.0;
}

bool WedgeBound::FindAngularBounds() {
  // Start with position directly above center point.
  AngularCoordinate start_ang;

  switch (sphere_) {
  case AngularCoordinate::Survey:
    start_ang.SetSurveyCoordinates(center_point_.Lambda()+radius_,
				   center_point_.Eta());
    break;
  case AngularCoordinate::Equatorial:
    start_ang.SetEquatorialCoordinates(center_point_.RA(),
				       center_point_.DEC()+radius_);
    break;
  case AngularCoordinate::Galactic:
    start_ang.SetGalacticCoordinates(center_point_.GalLon(),
				     center_point_.GalLat()+radius_);
    break;
  }

  double lammin = center_point_.Lambda(), lammax = center_point_.Lambda();
  double etamin = center_point_.Eta(), etamax = center_point_.Eta();

  AngularCoordinate rotate_ang;
  start_ang.Rotate(center_point_, position_angle_min_, rotate_ang, sphere_);
  if (DoubleLT(rotate_ang.Lambda(), lammin)) lammin = rotate_ang.Lambda();
  if (DoubleGT(rotate_ang.Lambda(), lammax)) lammax = rotate_ang.Lambda();
  if (DoubleLT(rotate_ang.Eta(), etamin)) etamin = rotate_ang.Eta();
  if (DoubleGT(rotate_ang.Eta(), etamax)) etamax = rotate_ang.Eta();

  start_ang.Rotate(center_point_, position_angle_max_, rotate_ang, sphere_);
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
  double circle_area = (1.0 - cos(radius_*DegToRad))*2.0*Pi*StradToDeg;
  SetArea(circle_area*(position_angle_max_ - position_angle_min_)/360.0);
  return true;
}

bool WedgeBound::CheckPoint(AngularCoordinate& ang) {
  bool within_bound = false;

  if (DoubleGE(center_point_.DotProduct(ang), costhetamin_)) {
    double cosphi = center_point_.CosPositionAngle(ang, sphere_);
    double sinphi = center_point_.SinPositionAngle(ang, sphere_);
    if (DoubleGE(cosphi, cosphi_min_) && DoubleLE(cosphi, cosphi_max_) &&
	DoubleGE(sinphi, sinphi_min_) && DoubleLE(sinphi, sinphi_max_))
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
    std::cout << "Stomp::PolygonBound::FindArea - " <<
      "Polygon area is over half the sphere.  This is bad.\n";
    return false;
  }

  if (tmp_area < 0.0) {
    std::cout << "Stomp::PolygonBound::FindArea - " <<
      "Polygon area is negative.  This is bad.\n";
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

  if (bound_area < Area() && Area() > 2.0*Pi*StradToDeg) {
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
