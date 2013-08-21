// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the ScalarMap class.  Unlike Maps, the primary
// goal here is to encode a scalar field over some area of the sky.  As such,
// we sacrifice some degree of precision in describing the exact area of the
// field and we use a uniform sampling of the field across the area in
// question.  This makes the class ideal for calculating angular correlation
// functions on the encoded field.


#ifndef SCALAR_UNION_H_
#define SCALAR_UNION_H_

#include <stdint.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <algorithm>
#include "angular_bin-inl.h"
#include "angular_correlation.h"
#include "core.h"
#include "point.h"
#include "scalar_pixel.h"
#include "bound_interface.h"
#include "circle_bound.h"

namespace s2omp {

class scalar_union: public bound_interface {
public:
  enum ScalarType {
    SCALAR_FIELD, DENSITY_FIELD, SAMPLED_FIELD
  };

  explicit scalar_union(ScalarType scalar_type);
  ~scalar_union();

  static scalar_union* from_bound(const bound_interface& bound, int level,
      ScalarType scalar_type);
  static scalar_union* from_scalar_pixels(const scalar_vector& pixels,
      ScalarType scalar_type);
  static scalar_union* from_scalar_union(const scalar_union& s, int level);

  bool init(const bound_interface& bound, int level, ScalarType t);
  bool init(const scalar_vector& pixels, ScalarType t);
  bool init(const scalar_union& u, int level);

  bool add_point(point& p, double intensity);
  bool add_point(point& p);

  scalar_pixel resample(const pixel& pix) const;

  double find_intensity(const pixel& pix) const;
  double find_density(const pixel& pix) const;
  double find_point_density(const pixel& pix) const;

  double find_local_area(const annulus_bound& bound) const;
  double find_local_intensity(const annulus_bound& bound) const;
  double find_local_density(const annulus_bound& bound) const;
  double find_local_point_density(const annulus_bound& bound) const;

  void calculate_mean_intensity();
  void convert_to_over_density();
  void convert_from_over_density();

  bool auto_correlate(angular_bin* theta);
  bool auto_correlate(angular_correlation* wtheta);

  bool auto_correlate_with_regions(const region_map& regions,
      angular_bin* theta);
  bool auto_correlate_with_regions(const region_map& regions,
      angular_correlation* wtheta);

  bool cross_correlate(scalar_union& s, angular_bin* theta);
  bool cross_correlate(scalar_union& s, angular_correlation* wtheta);

  bool cross_correlate_with_regions(scalar_union& s, const region_map& regions,
      angular_bin* theta);
  bool cross_correlate_with_regions(scalar_union& s, const region_map& regions,
      angular_correlation* wtheta);

  // Basic accessors for scalar_union parameters.
  inline int level() const {
    return level_;
  }
  inline ScalarType type() const {
    return type_;
  }
  inline double intensity() const {
    return total_intensity_;
  }
  inline int n_points() const {
    return total_points_;
  }
  inline double density() const {
    return total_points_ == 0 ? total_intensity_ : total_intensity_
        / total_points_;
  }
  inline double point_density() const {
    return total_points_ / area_;
  }

  inline double mean_intensity() const {
    return mean_intensity_;
  }
  inline bool is_over_density() const {
    return converted_to_overdensity_;
  }
  inline scalar_const_iterator begin() const {
    return pixels_.begin();
  }
  inline scalar_const_iterator end() const {
    return pixels_.end();
  }

  // API from pixelized_bound_interface.h
  virtual bool is_empty() const {
    return pixels_.empty();
  }

  virtual long size() const {
    return pixels_.size();
  }

  virtual void clear();

  virtual double area() const {
    return area_;
  }

  virtual bool contains(const point& p) const;
  virtual bool contains(const pixel& pix) const;

  virtual double contained_area(const pixel& pix) const;
  virtual bool may_intersect(const pixel& pix) const;

  virtual circle_bound get_bound() const;
  virtual point get_center() const;

  virtual void get_covering(pixel_vector* pixels) const;
  virtual void get_covering(const long max_pixels, pixel_vector* pixels) const;
  virtual void get_covering(double fractional_area_tolerance,
      pixel_vector* pixels) const;
  virtual void get_interior_covering(int max_level, pixel_vector* pixels) const;
  virtual void get_simple_covering(int level, pixel_vector* pixels) const;
  virtual void get_center_covering(int level, pixel_vector* pixels) const;

private:
  scalar_union();
  void initialize_bound();
  bool correlate_unions(scalar_union& s, const region_map& regions,
      bool autocorrelate, angular_bin* theta);

  scalar_vector pixels_;
  double area_, mean_intensity_, unmasked_fraction_minimum_, total_intensity_;
  ScalarType type_;
  int level_;
  uint32_t total_points_;
  circle_bound bound_;
  bool initialized_, initialized_bound_;
  bool converted_to_overdensity_, calculated_mean_intensity_;
};

} // end namespace s2omp

#endif /* SCALAR_UNION_H_ */
