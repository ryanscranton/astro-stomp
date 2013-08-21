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

#include "core.h"

namespace s2omp {

class bound_interface {
public:
  virtual ~bound_interface();

  // All bound objects should implement versions of these basic methods.  We
  // don't expect that all of them will necessarily make sense for all
  // derived classes (e.g., what is the size of a circle_bound?), but enough
  // should to justify each method.  We start with the basics: is a bound
  // empty and, if not, what's its area.
  virtual bool is_empty() const;
  virtual long size() const;
  virtual void clear();
  virtual double area() const;

  // As bound-derived objects represent regions on the sky, we should be able
  // to do some common tasks with them.  In particular, we should be able to
  // create a low-resolution covering of the bound that gives an approximation
  // of the area.  Likewise, we should be able to split the area covered by
  // the bound into roughly equal-area pieces.  Finally, we should be able to
  // generate random points over the bound's extent.

  // For these tasks, we need a few common methods whose implementations will
  // vary based on the details of the different bound objects.  First, we need
  // to be able to know if a pixel or point is contained within the bound.
  virtual bool contains(const point& p) const;
  virtual bool contains(const pixel& pix) const;

  // For regionation, we need to know how much of a given pixel is contained
  // within the bound.
  virtual double contained_area(const pixel& pix) const;

  // Since generating a covering is somewhat inexact, we only need to know if
  // a pixel might intersect our bound, which is generally something we can
  // calculate faster than containment.
  virtual bool may_intersect(const pixel& pix) const;

  // To facilitate generating a covering, we need a starting point.  This is
  // generally going to be something like a centroid, but should be guaranteed
  // to be contained within the bound if at all possible.
  virtual point get_center() const;

  // For the purposes of generating random points, we need to be able
  // to create a circle_bound that encompasses the bound.
  virtual circle_bound get_bound() const;

  // For all derived instances of bound_interface, we implement a generic
  // method for generating random points within that object's bound.  These
  // should be over-ridden if the details of that object's implementation make
  // calling get_bound() expensive.
  virtual point get_random_point();
  virtual void get_random_points(long n_points, point_vector* points);

  virtual void get_covering(pixel_vector* pixels) const;
  virtual void get_covering(
      long max_pixels, pixel_vector* pixels) const;
  virtual void get_covering(
      double fractional_area_tolerance, pixel_vector* pixels) const;
  virtual void get_interior_covering(int max_level, pixel_vector* pixels) const;
  virtual void get_simple_covering(int level, pixel_vector* pixels) const;
  virtual void get_center_covering(int level, pixel_vector* pixels) const;
};

} // end namespace s2omp

#endif /* BOUND_INTERFACE_H_ */
