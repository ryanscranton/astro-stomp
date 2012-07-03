// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the abstract BaseMap class that serves as the
// basis for all of the *Map objects.  BaseMap sets out the basic functionality
// that all of the *Map classes need to describe a given region on the sky and
// do some basic internal maintenance.  Additionally, BaseMap provides a
// common set of methods for dividing that area up into nearly equal-area,
// similarly-shaped regions.  This functionality is the basis for calculating
// jack-knife errors for our various statistical analyses.

#ifndef BASE_BOUND_H_
#define BASE_BOUND_H_

#include <stdint.h>
#include <vector>
#include <map>
#include "stomp_core.h"
#include "stomp_geometry.h"
#include "stomp_pixel.h"

namespace s2omp {

class point;  // class declaration in stomp_angular_coordinate.h
class pixel;              // class declaration in stomp_pixel.h
class circle_bound;
class base_bound;

class base_bound {
public:
  virtual bool is_empty();
  virtual long size();
  virtual void clear();
  virtual void area();

  virtual bool contains(point& p);
  virtual bool contains(pixel& pix);
  virtual bool contains(base_bound& b);

  virtual double contained_area(pixel& pix);
  virtual bool may_intersect(pixel& pix);

  virtual void covering(pixel_vector& pixels);
  virtual void covering(int max_pixels, pixel_vector& pixels);
  virtual void covering(double fractional_area_tolerance, pixel_vector& pixels);
  virtual void simple_covering(int level, pixel_vector& pixels);

  virtual circle_bound get_bound();

  virtual point get_random_point();
  virtual void get_random_points(long n_points, pixel_vector& points);

private:
  base_bound();
  virtual ~base_bound();
};


} // end namespace s2omp
#endif /* BASE_BOUND_H_ */
