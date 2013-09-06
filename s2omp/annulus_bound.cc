/*
 * annulus_bound.cc
 *
 *  Created on: Aug 13, 2012
 *      Author: cbmorrison
 */

#include "angular_bin-inl.h"
#include "annulus_bound.h"
#include "bound_interface.h"
#include "circle_bound.h"
#include "pixel.h"

namespace s2omp {

annulus_bound::annulus_bound() {
  axis_ = point(1.0, 0.0, 0.0, 1.0);
  inner_height_ = -1.0;
  outer_height_ = -1.0;
  update_caps();
}

annulus_bound::annulus_bound(const point& axis, double inner_height,
    double outer_height) {
  axis_ = axis;
  inner_height_ = inner_height;
  outer_height_ = outer_height;
  update_caps();
}

annulus_bound::~annulus_bound() {
}

annulus_bound* annulus_bound::from_radii(const point& axis,
    double inner_radius_degrees, double outer_radius_degrees) {
  return new annulus_bound(axis, circle_bound::get_height_for_angle(
      inner_radius_degrees), circle_bound::get_height_for_angle(
      outer_radius_degrees));
}

annulus_bound* annulus_bound::from_heights(const point& axis,
    double inner_height, double outer_height) {
  return new annulus_bound(axis, inner_height, outer_height);
}

annulus_bound* annulus_bound::from_angular_bin(const point& axis,
    const angular_bin& bin) {
  return from_heights(axis,
      circle_bound::get_height_for_angle(bin.theta_min()),
      circle_bound::get_height_for_angle(bin.theta_max()));
}

bool annulus_bound::is_empty() const {
  return !is_valid();
}

long annulus_bound::size() const {
  return is_empty() ? 0 : 1;
}

void annulus_bound::clear() {
  axis_ = point();
  inner_height_ = -1.0;
  outer_height_ = -1.0;
}

double annulus_bound::area() const {
  return is_empty() ? 0.0 : 2.0 * PI * (outer_height_ - inner_height_)
      * STRAD_TO_DEG2;
}

bool annulus_bound::contains(const point& p) const {
  return outer_contains(p) && !inner_contains(p);
}

bool annulus_bound::contains(const pixel& pix) const {
  // This one works for the complement case
  return outer_contains(pix) && !inner_may_intersect(pix);
}

double annulus_bound::contained_area(const pixel& pix) const {
  if (contains(pix)) {
    return pix.exact_area();
  }

  if (!may_intersect(pix)) {
    return 0.0;
  }

  // TODO(scranton): If the pixel is considerably larger than our bound, then
  // we can't efficiently do sub-sampling to figure out the overlapping area.
  // In that case, since we know that, at this point, the bound intersects
  // the pixel, then the most likely case is that the bound is entirely
  // contained by pixel.  This is a bit of a kludge, but lacking the ability to
  // do this calculation analytically, we just return the area of the bound.
  if (pix.exact_area() > 10.0 * area()) {
    return area();
  }

  // TODO(scranton): This is something that could probably be done analytically,
  // but for now, we'll do something simpler.
  double total_area = 0.0;
  int sampling_level = min(pix.level() + 5, MAX_LEVEL);
  pixel child_end = pix.child_end(sampling_level);
  for (pixel child_pix = pix.child_begin(sampling_level); child_pix
      != child_end; child_pix = child_pix.next()) {
    total_area += contains(child_pix) ? child_pix.exact_area() : 0.0;
  }

  return total_area;
}

bool annulus_bound::may_intersect(const pixel& pix) const {
  return outer_may_intersect(pix) && !inner_contains(pix);
}

} // end namespace s2omp
