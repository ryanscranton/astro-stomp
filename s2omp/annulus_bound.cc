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
  return new annulus_bound(
      axis, circle_bound::get_height_for_angle(inner_radius_degrees),
      circle_bound::get_height_for_angle(outer_radius_degrees));
}

annulus_bound* annulus_bound::from_heights(const point& axis,
    double inner_height, double outer_height) {
  return new annulus_bound(axis, inner_height, outer_height);
}

annulus_bound* annulus_bound::from_angular_bin(const point& axis,
    const angular_bin& bin) {
  return from_heights(
      axis, circle_bound::get_height_for_angle(bin.theta_min()),
      circle_bound::get_height_for_angle(bin.theta_max()));
}

bool annulus_bound::is_empty() const {
  return outer_bound_->is_empty();
}

long annulus_bound::size() const {
  return is_empty() ? 0 : 2;
}

void annulus_bound::clear() {
  delete outer_bound_;
  delete inner_bound_;
}

double annulus_bound::area() const {
  return is_empty() ? 0.0 : outer_bound_->area() - inner_bound_->area();
}

bool annulus_bound::contains(const point& p) const {
  return (outer_bound_->contains(p) && !inner_bound_->contains(p));
}

bool annulus_bound::contains(const pixel& pix) const {
  // This one works for the complement case
  return (outer_bound_->contains(pix) && !inner_bound_->may_intersect(pix));
}

double annulus_bound::contained_area(const pixel& pix) const {
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

bool annulus_bound::may_intersect(const pixel& pix) const {
  return (outer_bound_->may_intersect(pix) &&
          !inner_bound_->contains(pix));
}

bool annulus_bound::may_intersect_outer(const pixel& pix) const {
  return outer_bound_->may_intersect(pix);
}

} // end namespace s2omp
