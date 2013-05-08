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
#include "pixel.h";

namespace s2omp {

annulus_bound::annulus_bound(const point& axis, double inner_height,
                             double outer_height) {
  outer_bound_ = circle_bound::from_height(axis, outer_height);
  inner_bound_ = circle_bound::from_height(axis, inner_height);
}

annulus_bound::~annulus_bound() {
  delete outer_bound_;
  delete inner_bound_;
}

annulus_bound* annulus_bound::from_radii(
    const point& axis, double inner_radius_degrees,
    double outer_radius_degrees) {
  return new annulus_bound(axis, 1.0 - cos(inner_radius_degrees * DEG_TO_RAD),
                           1.0 - cos(outer_radius_degrees * DEG_TO_RAD));
}

annulus_bound* annulus_bound::from_heights(const point& axis,
    double inner_height, double outer_height) {
  return new annulus_bound(axis, inner_height, outer_height);
}

annulus_bound* annulus_bound::from_angular_bin(const point& axis,
    const angular_bin& bin) {
  return from_heights(
      axis, 1.0 - bin.cos_theta_min(), 1.0 - bin.cos_theta_max());
}

bool annulus_bound::is_empty() {
  return outer_bound_->is_empty();
}

long annulus_bound::size() {
  if (!is_empty())
    return 2;
  return 0;
}

void annulus_bound::clear() {
  delete outer_bound_;
  delete inner_bound_;
}

double annulus_bound::area() {
  if (!is_empty())
    return outer_bound_->area() - inner_bound_->area();
  return 0.0;
}

bool annulus_bound::contains(const point& p) {
  return (outer_bound_->contains(p) && !inner_bound_->contains(p));
}

bool annulus_bound::contains(const pixel& pix) {
  // This one works for the complement case
  return (outer_bound_->contains(pix) && !inner_bound_->may_intersect(pix));
}

double annulus_bound::contained_area(const pixel& pix) {
 double area = 0.0;
 if (contains(pix)) {
   return pix.exact_area();
 } else if (may_intersect(pix)) {
   for (pixel c = pix.child_begin(); c != pix.child_end(); c = c.next()) {
     area += contained_area(c);
   }
 }
 return area;
}

bool annulus_bound::may_intersect(const pixel& pix) {
  return (outer_bound_->may_intersect(pix) &&
          !inner_bound_->contains(pix));
}

bool annulus_bound::may_intersect_outer(const pixel& pix) {
  return outer_bound_->may_intersect(pix);
}

} // end namespace s2omp
