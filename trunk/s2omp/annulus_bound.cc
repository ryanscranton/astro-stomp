/*
 * annulus_bound.cc
 *
 *  Created on: Aug 13, 2012
 *      Author: cbmorrison
 */




#include "annulus_bound.h"

namespace s2omp {

annulus_bound::annulus_bound() {
  inner_height_ = -1;
  outer_height_ = -1;
}

annulus_bound::annulus_bound(const point& axis, double inner_height,
    double outer_height) {
  axis_ = axis;
  inner_height_ = inner_height;
  outer_height_ = outer_height;
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

annulus_bound* annulus_bound::from_angular_bin(const point& axis,
    const angular_bin& bin) {
  return from_radii(axis, bin.theta_min(), bin.theta_max());
}

bool annulus_bound::is_empty() {
  if (inner_height_ < 0 || outer_height_ < 0)
    return true;
  return false;
}

long annulus_bound::size() {
  if (!is_empty())
    return 2;
  return 0;
}

double annulus_bound::area() {
  if (!is_empty_())
    return 2.0 * PI * (outer_height_ - inner_height_)*STRAD_TO_DEG2;
  return 0.0;
}

bool annulus_bound::contains(const point& p) {
  double p_height = 1.0 - axis_.dot(p);
  if (outer_height_ >= p_height && p_height >= inner_height_)
    return true;
  return false;
}

bool annulus_bound::contains(const pixel& pix) {
  // This one works for the complement case
  if (contains(pix, outer_height_) && !may_intersect(pix, inner_height_))
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
  if (may_intersect(pix, outer_height_) && !contains(pix, inner_height_))
    return true;
  return false;
}

bool annulus_bound::may_intersect_outer(const pixel& pix) {
  return may_intersect(pix, outer_height_);
}

circle_bound* annulus_bound::get_bound() {
  return circle_bound.from_height(axis_, outer_height_);
}

bool annulus_bound::contains(const pixel& pix, double height) {
  for (int k = 0; k < 4; ++k) {
    if (!contains(pix.vertex(k))) return false;
  }
  circle_bound* comp = complement(height);
  bool ans = !comp->may_intersect(pix);
  delete comp;
  return ans;
}

bool annulus_bound::may_intersect(const pixel& pix, double height) {
  point_vector vertices;
  vertices.reserve(4);
  for (int k = 0; k < 4; ++k) {
    vertices.push_back(pix.vertex(4));
    if (contains(vertices.back())) return true;
  }
  return intersects(pix, vertices, height);
}

bool annulus_bound::intersects(const pixel& pix, point_vector vertices,
    double height) {
  // Much of this code is lifted from S2::S2Cap.Intersects

  // If the cap that we are considering is a hemisphere or larger, then since
  // there are no cells that would intersect this cell that would not already
  // have a vertex contained.
  if (height >= 1) return false;

  // Check for empty caps before checking if the axis is contained by the
  // pixel.
  if (is_empty()) return false;

  // If there are no vertices contained in the bound, but the axis is contained
  // then this bound intersects the pixel.
  if (pix.contains(axis_)) return true;

  // Since the vertices aren't contained and the axis isn't contained then
  // the only points left to test are those long the interior of an edge.

  double sin2_angle = height * (2 - height);
  for (int k = 0; k < 4; ++k) {
    point edge = pix.edge(k);
    double dot = axis_.dot(edge);
    if (dot > 0) {
      // If the dot product is positive then the the current edge is not the
      // edge of closest approach and, since no vertices are contained, we
      // don't need to consider it.
      continue;
    }

    if (dot * dot > sin2_angle){
      // If this is the case then the closest point on the edge to the axis
      // is outside that of the circle_bound that defines this region. We then
      // return false.
      return false;
    }

    // Since we've gone through the above tests we know that the edge passes
    // through the bound at some point along the great cirlce. We now need to
    // test if this point is between the two vertices of the edge.
    point dir = edge.cross(axis_);
    if (dir.dot(vertices[k]) < 0 && dir.dot(vertices[(k+1)&3]) > 0)
      return true;
  }
  return false;
}

circle_bound* annulus_bound::complement(double height) {
  circle_bound* bound = circle_bound.from_height(-axis_, 2.0 - height);
  return bound;
}

} // end namespace s2omp
