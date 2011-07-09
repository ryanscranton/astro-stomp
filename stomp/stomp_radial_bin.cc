// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//

#include "stomp_core.h"
#include "stomp_angular_bin.h"
#include "stomp_radial_bin.h"
#include "stomp_util.h"

namespace Stomp {

RadialBin::RadialBin() {
  r_min_ = r_max_ = redshift_ = 0.0;
}

RadialBin::RadialBin(double r_min, double r_max, double redshift) {
  r_min_ = r_min;
  r_max_ = r_max;
  redshift_ = redshift;
  SetThetaMin(Cosmology::ProjectedAngle(redshift_, r_min_));
  SetThetaMax(Cosmology::ProjectedAngle(redshift_, r_max_));
}

RadialBin::RadialBin(double r_min, double r_max, double redshift,
                     int16_t n_regions) {
  r_min_ = r_min;
  r_max_ = r_max;
  redshift_ = redshift;
  SetThetaMin(Cosmology::ProjectedAngle(redshift_, r_min_));
  SetThetaMax(Cosmology::ProjectedAngle(redshift_, r_max_));
  ClearRegions();
  if (n_regions > 0) InitializeRegions(n_regions);
}

RadialBin::~RadialBin() {
  r_min_ = r_max_ = redshift_ = 0.0;
  ClearRegions();
}

void RadialBin::SetRadius(double r) {
  r_ = r;
  SetTheta(Cosmology::ProjectedAngle(redshift_, r_));
}

void RadialBin::SetRadiusMin(double r_min) {
  r_min_ = r_min;
  SetThetaMin(Cosmology::ProjectedAngle(redshift_, r_min_));
}

void RadialBin::SetRadiusMax(double r_max) {
  r_max_ = r_max;
  SetThetaMax(Cosmology::ProjectedAngle(redshift_, r_max_));
}

void RadialBin::SetRedshift(double z) {
  redshift_ = z;
  SetThetaMin(Cosmology::ProjectedAngle(redshift_, r_min_));
  SetThetaMax(Cosmology::ProjectedAngle(redshift_, r_max_));
}

bool RadialBin::WithinRadialBounds(double r) {
  return (DoubleGE(r, r_min_) &&
	  DoubleLE(r, r_max_) ? true : false);
}

double RadialBin::Radius() {
  return r_;
}

double RadialBin::RadiusMin() {
  return r_min_;
}

double RadialBin::RadiusMax() {
  return r_max_;
}

double RadialBin::Redshift() {
  return redshift_;
}

bool RadialBin::RadialOrder(RadialBin r_a, RadialBin r_b) {
  return (r_a.ThetaMin() < r_b.ThetaMin() ? true : false);
}

bool RadialBin::ReverseResolutionOrder(RadialBin r_a,
                                       RadialBin r_b) {
  return (r_b.Resolution() < r_a.Resolution() ? true : false);
}

} // end namespace Stomp
