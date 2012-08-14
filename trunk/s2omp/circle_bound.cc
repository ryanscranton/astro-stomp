/*
 * circle_bound.cc
 *
 *  Created on: Jul 18, 2012
 *      Author: cbmorrison
 */




#include "circle_bound.h"

namespace s2omp {

circle_bound::circle_bound() {
  cap_ = S2::S2Cap();
}

circle_bound::circle_bound(const point& axis, double height) {
  cap_ = S2::S2Cap.FromAxisHeight(axis.s2point(), height);
}

static circle_bound* cirlce_bound::from_angular_bin(const point& axis,
      const angular_bin& bin) {
  circle_bound* circle = new circle_bound(axis,
      1.0 - cos(bin.theta_max()*DEG_TO_RAD));
  return circle;
}

circle_bound* circle_bound::from_radius(point& axis, double radius_degrees) {
  double height = 1.0 - cos(radius_degrees*DEG_TO_RAD);
  circle_bound* circle = new circle_bound(axis, height);
  return circle;
}

cirlce_bound* circle_bound::from_height(point& axis, double height) {
  circle_bound* circle = new circle_bound();
  circle_bound->cap_ = S2::S2Cap.FromAxisHeight(axis.s2point(),height);
  return circle;
}

point circle_bound::axis() {
  return point.from_s2point(cap_.axis(), 1.0);
}

virtual bool circle_bound::is_empty() {
  return cap_.is_empty();
}

virtual long circle_bound::size() {
  if (is_empty())
    return 0;
  return 1;
}

virtual double circle_bound::area() {
  return cap_.area()*STRAD_TO_DEG2;
}

virtual circle_bound circle_bound::get_bound() {
  return *this;
}

bool circle_bound::contains(const point& p) {
  return cap_.Contains(p.s2point());
}

// TODO(cbmorrison) the intermediate S2Cell creation may be expensive and slow.
// If needed create an s2omp class instead of wrapping cap_.contains.
bool circle_bound::contains(const pixel& pix) {
  return cap_.Contains(S2::S2Cell(pix.id()));
}

double circle_bound::contained_area(const pixel& pix) {
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

// TODO(cbmorrison) the intermediate S2Cell creation may be expensive and slow.
// If needed create an s2omp class instead of wrapping cap_.may_intersect.
bool circle_bound::may_intersect(const pixel& pix) {
  return cap_.MayIntersect(S2::S2Cell(pix.id()));
}

virtual point circle_bound::get_random_point() {
  if (!initialized_random_) {
    initialize_random();
  }
  // To generate a random point on the sphere we first generate the point as
  // if this cap's axis is the z axis. We then rotate the point to the true
  // position on the sphere.

  // Generate our random values. We first generate a point uniform in
  // cos(theta) from the lowest point on the cap to the highest. Next we
  // generate a point uniform in phi.
  double z = mtrand.rand(cap_.height()) + 1 - cap_.height();
  double phi = mtrand.rand(2.0*PI);

  // To turn our angles into x,y,z coordinates we need to now the sin as well.
  double sintheta = sin(acos(z));

  // Now that we have the random point we need to rotate it to the correct
  // position on the sphere. We do this by rotating the point around the normal
  // to the great circle defined by the z axis and the cap axis. We rotate by
  // the angle the cap axis makes with the z axis.
  point p = point(sintheta*cos(phi), sintheta*sin(phi), z);
  p.rotate_about(great_circle_norm_, rotate_);

  return p;
}

void void circle_bound::get_random_points(long n_points,
    point_vector* points) {
  if (!points->empty())
    points->clear();
  points->reserve(n_points);

  for (long i = 0; i < n_points; ++i) {
    points->push_back(get_random_point());
  }
}

point circle_bound::get_weighted_random_point(const point_vector& points) {
  point p = get_random_point();
  p.set_weight(points[mtrand.randInt(points.size())].weight());
  return p;
}

void circle_bound::get_weighted_random_points(long n_points,
    point_vector* points, const point_vector& input_points) {
  if (!points->empty()) points->clear();

  for (long i = 0; i < n_points; ++i) {
    point p = get_random_point();
    p.set_weight(mtrand.randInt(input_points.size()));
    points->push_back(p);
  }
}

void circle_bound::initialize_random() {
  mtrand.seed();
  rotate_ = acos(cap_.axis()[2]);
  great_circle_norm_ = point(-cap_.axis[1], cap_.axis[0], 0.0, 1.0);
}

} //end namspace s2omp
