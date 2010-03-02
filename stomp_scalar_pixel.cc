// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains a variant on the Pixel class.  The goal here is to
// encode some manner of scalar field (galaxy density on the sky, CMB
// temperature, etc.) where we need both the value and an associated noise on
// the value.  Likewise, as will be seen in the ScalarMap class, these pixels
// will generally be used to sample this scalar field uniformly over some
// region.  This is in contrast with the Map object where the goal is to
// accurately describe a region's geometry using pixels of various sizes.

#include "stomp_core.h"
#include "stomp_scalar_pixel.h"

namespace Stomp {

ScalarPixel::ScalarPixel() {
  intensity_ = 0.0;
  n_point_ = 0;
}

ScalarPixel::ScalarPixel(const uint16_t input_resolution,
			 const uint32_t input_pixnum,
			 const double input_weight,
			 const double input_intensity,
			 const uint32_t n_point) {
  SetResolution(input_resolution);

  uint32_t tmp_y = input_pixnum/(Nx0*Resolution());
  uint32_t tmp_x = input_pixnum - Nx0*Resolution()*tmp_y;

  SetPixnumFromXY(tmp_x, tmp_y);
  SetWeight(input_weight);
  intensity_ = input_intensity;
  n_point_ = n_point;
}

ScalarPixel::ScalarPixel(const uint16_t input_resolution,
			 const uint32_t input_hpixnum,
			 const uint32_t input_superpixnum,
			 const double input_weight,
			 const double input_intensity,
			 const uint32_t n_point) {
  SetResolution(input_resolution);

  uint32_t tmp_x, tmp_y;
  Pixel::HPix2XY(input_resolution, input_hpixnum, input_superpixnum,
		 tmp_x, tmp_y);

  SetPixnumFromXY(tmp_x, tmp_y);
  SetWeight(input_weight);
  intensity_ = input_intensity;
  n_point_ = n_point;
}

ScalarPixel::ScalarPixel(const uint32_t input_x,
			 const uint32_t input_y,
			 const uint16_t input_resolution,
			 const double input_weight,
			 const double input_intensity,
			 const uint32_t n_point) {
  SetResolution(input_resolution);
  SetPixnumFromXY(input_x, input_y);
  SetWeight(input_weight);
  intensity_ = input_intensity;
  n_point_ = n_point;
}

ScalarPixel::ScalarPixel(AngularCoordinate& ang,
			 const uint16_t input_resolution,
			 const double input_weight,
			 const double input_intensity,
			 const uint32_t n_point) {

  SetResolution(input_resolution);
  SetPixnumFromAng(ang);
  SetWeight(input_weight);
  intensity_ = input_intensity;
  n_point_ = n_point;
}

ScalarPixel::~ScalarPixel() {
  intensity_ = 0.0;
  n_point_ = 0;
}

void ScalarPixel::SetIntensity(const double intensity) {
  intensity_ = intensity;
}

void ScalarPixel::SetNPoints(const uint32_t n_point) {
  n_point_ = n_point;
}

double ScalarPixel::Intensity() {
  return intensity_;
}

uint32_t ScalarPixel::NPoints() {
  return n_point_;
}

double ScalarPixel::MeanIntensity() {
  return (n_point_ == 0 ? intensity_ : intensity_/n_point_);
}

void ScalarPixel::AddToIntensity(const double intensity,
				 const uint32_t n_point) {
  intensity_ += intensity;
  n_point_ += n_point;
}

void ScalarPixel::ScaleIntensity(double scale_factor) {
  intensity_ *= scale_factor;
}

void ScalarPixel::NormalizeIntensity() {
  intensity_ /= n_point_;
  n_point_ = 0;
}

void ScalarPixel::ConvertToOverDensity(double expected_intensity) {
  intensity_ -= expected_intensity;
  is_overdensity_ = true;
}

void ScalarPixel::ConvertToFractionalOverDensity(double expected_intensity) {
  intensity_ =
    (intensity_ - expected_intensity*Weight()*Area())/
    (expected_intensity*Weight()*Area());
  is_overdensity_ = true;
}

void ScalarPixel::ConvertFromOverDensity(double expected_intensity) {
  intensity_ += expected_intensity;
  is_overdensity_ = false;
}

void ScalarPixel::ConvertFromFractionalOverDensity(double expected_intensity) {
  double norm_intensity = expected_intensity*Weight()*Area();
  intensity_ = intensity_*norm_intensity + norm_intensity;
  is_overdensity_ = false;
}

bool ScalarPixel::IsOverDensity() {
  return is_overdensity_;
}

} // end namespace Stomp
