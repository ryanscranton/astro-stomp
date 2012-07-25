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

#ifndef BOUND_INTERFACE_H_
#define BOUND_INTERFACE_H_

namespace s2omp {

class bound_interface {
public:
  virtual ~bound_interface();

  virtual bool is_empty();
  virtual long size();
  virtual void clear();
  virtual void area();

  virtual bool contains(const point& p);
  virtual bool contains(const pixel& pix);
  virtual bool contains(const bound& b);

  virtual double contained_area(const pixel& pix);
  virtual bool may_intersect(const pixel& pix);

  virtual void covering(pixel_vector* pixels);
  virtual void covering(int max_pixels, pixel_vector* pixels);
  virtual void covering(double fractional_area_tolerance, pixel_vector* pixels);
  virtual void simple_covering(int level, pixel_vector* pixels);

  virtual circle_bound get_bound();

  virtual point get_random_point();
  virtual void get_random_points(long n_points, pixel_vector* points);

private:
  bound_interface();
};

} // end namespace s2omp

#endif /* _H_ */
