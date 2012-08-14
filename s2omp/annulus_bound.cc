/*
 * annulus_bound.cc
 *
 *  Created on: Aug 13, 2012
 *      Author: cbmorrison
 */




#include "annulus_bound.h"

namespace s2omp {

annulus_bound::annulus_bound() {
  inner_cap_ = S2::S2Cap();
  outer_cap_ = S2::S2Cap();
}

annulus_bound::annulus_bound(const point& axis, double inner_height,
    double outer_costheta) {
  inner_cap_ = S2::S2Cap::FromAxisHeight(axis.s2point(), inner_height);
  outer_cap_ = S2::S2Cap::FromAxisHeight(axis.s2point(), outer_height);
}

annulus_bound* annulus_bound::from_radii(const point& axis,
    double inner_radius_degrees, double outer_radius_degrees) {
  annulus_bound* bound = new annulus_bound(axis,
      1.0 - cos(inner_radius_degrees*DEG_TO_RAD),
      1.0 - cos(outer_radius_degrees*DEG_TO_RAD));
  return bound;
}

annulus_bound* annulus_bound::from_heights(const point& axis,
    double inner_height, double outer_height) {
  annulus_bound* bound = new annulus_bound(axis, inner_height, outer_height);
  return bound;
}

bool annulus_bound::is_empty() {
  if (inner_cap_.is_empty() || outer_cap_.is_empty())
    return true;
  return false;
}

long annulus_bound::size() {
  if (!is_empty())
    return 2;
  return 0;
}

double annulus_bound::area() {
  return (outer_cap_.area() - inner_cap_.area())*STRAD_TO_DEG2;
}

bool annulus_bound::contains(const point& p) {
  if (outer_cap_.contains(p.s2point()) && !inner_cap_.contains(p.s2point()))
    return true;
  return false;
}

bool annulus_bound::contains(const pixel& pix) {
  if (outer_cap_.contains(S2::S2Cell(pix.id())) &&
      !inner_cap_.MayIntersect(S2::S2Cell(pix.id())))
    return true;
  return false;
}

double annulus_bound::contained_area(const pixel& pix) {
 double area = 0.0;
 if (contains(pix)) {
   return pix.exact_area();
 } else if (may_intersect(pix)) {
   for (pixel_iterator iter = pix.child_begin();
       iter != pix.child_end(); ++iter) {
     area += contained_area(*iter);
   }
 }
 return area;
}

bool annulus_bound::may_intersect(const pixel& pix) {
  if (outer_cap_.MayIntersect(S2::S2Cell(pix.id())) &&
      !inter_cap_.Contains(S2::S2Cell(pix.id())))
    return true;
  return false;
}

circle_bound* annulus_bound::get_bound() {
  return cirlce_bound.from_height(point(outer_cap_.axis()),
      outer_cap_.height());
}



} // end namespace s2omp
