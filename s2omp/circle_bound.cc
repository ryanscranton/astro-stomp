/*
 * circle_bound.cc
 *
 *  Created on: Jul 18, 2012
 *      Author: cbmorrison
 */

#include <s2/s2cell.h>

#include "angular_bin-inl.h"
#include "circle_bound.h"
#include "pixel.h"
#include "point.h"

namespace s2omp {

circle_bound::circle_bound() {
  axis_ = point();
  height_ = -1.0;
}

circle_bound::circle_bound(const point& axis, double height) {
  axis_ = axis;
  height_ = height;
}

circle_bound::~circle_bound() {
}

double circle_bound::get_height_for_angle(double theta_degrees) {
  double theta_radians = theta_degrees * DEG_TO_RAD;
  if (theta_radians > PI) {
    return 2.0;
  }

  // Following S2's convention, we return 2 * sin^2(theta/2) instead of
  // 1.0 - cos(theta) since the former is more accurate than the latter for
  // theta close to 0 in radians.
  double d = sin(0.5 * theta_radians);
  return 2.0 * d * d;
}

circle_bound* circle_bound::from_angular_bin(const point& axis,
    const angular_bin& bin) {

  return new circle_bound(axis, get_height_for_angle(bin.theta_max()));
}

circle_bound* circle_bound::from_radius(const point& axis,
    double radius_degrees) {
  return new circle_bound(axis, get_height_for_angle(radius_degrees));
}

circle_bound* circle_bound::from_height(const point& axis, double height) {
  return new circle_bound(axis, height);
}

void circle_bound::add_point(const point& p) {
  if (is_empty()) {
    axis_ = p;
    height_ = 0.0;
  } else {
    // Calculate the chord distance between the input point and the central
    // axis.  Increase the height if this is larger than the current value.
    double dist2 = (axis_.unit_sphere_x() - p.unit_sphere_x())
        * (axis_.unit_sphere_x() - p.unit_sphere_x()) + (axis_.unit_sphere_y()
        - p.unit_sphere_y()) * (axis_.unit_sphere_y() - p.unit_sphere_y())
        + (axis_.unit_sphere_z() - p.unit_sphere_z()) * (axis_.unit_sphere_z()
            - p.unit_sphere_z());
    height_ = max(height_, FLOAT_ROUND_UP * 0.5 * dist2);
  }
}

void circle_bound::add_circle_bound(const circle_bound& bound) {
  if (is_empty()) {
    *this = bound;
  } else {
    // Calculate the maximum angular distance between the axis and the edge of
    // the input circle_bound.  Increase the height accordingly.
    double angle = axis_.angular_distance(bound.axis()) + bound.radius();
    height_ = max(height_, FLOAT_ROUND_UP * get_height_for_angle(angle));
  }
}

point circle_bound::create_random_point() {
  if (!initialized_random_) {
    initialize_random();
  }

  // To generate a random point on the sphere we first generate the point as
  // if this cap's axis is the z axis. We then rotate the point to the true
  // position on the sphere.

  // Generate our random values. We first generate a point uniform in
  // cos(theta) from the lowest point on the cap to the highest. Next we
  // generate a point uniform in phi.
  double z = mtrand_.rand(height_) + 1 - height_;
  double phi = mtrand_.rand(2.0 * PI);

  // To turn our angles into x,y,z coordinates we need to now the sin as well.
  double sintheta = sin(acos(z));

  // Now that we have the random point we need to rotate it to the correct
  // position on the sphere. We do this by rotating the point around the normal
  // to the great circle defined by the z axis and the cap axis. We rotate by
  // the angle the cap axis makes with the z axis.
  point p = point(sintheta * cos(phi), sintheta * sin(phi), z, 1.0);
  p.rotate_about(great_circle_norm_, rotate_);

  return p;
}

void circle_bound::get_weighted_random_points(long n_points,
    point_vector* points, const point_vector& input_points) {
  if (!points->empty())
    points->clear();

  for (long i = 0; i < n_points; ++i) {
    point p = create_random_point();
    p.set_weight(mtrand_.randInt(input_points.size()));
    points->push_back(p);
  }
}

bool circle_bound::is_empty() const {
  return height_ < 0;
}

long circle_bound::size() const {
  if (is_empty())
    return 0;
  return 1;
}

void circle_bound::clear() {
  height_ = 0.0;
  axis_ = point();
}

double circle_bound::area() const {
  return 2.0 * PI * height_ * STRAD_TO_DEG2;
}

bool circle_bound::contains(const point& p) const {
  double p_height = 1.0 - axis_.dot(p);
  return height_ >= p_height;
}

bool circle_bound::contains(const pixel& pix) const {
  return get_s2cap().Contains(pix.get_cell());
}

double circle_bound::contained_area(const pixel& pix) const {
  if (contains(pix)) {
    return pix.exact_area();
  }

  if (!may_intersect(pix)) {
    return 0.0;
  }

  // TODO(scranton): This is something that could probably be done analytically,
  // but for now, we'll do something simpler.
  double total_area = 0.0;
  int sampling_level = min(pix.level() + 5, MAX_LEVEL);
  pixel child_end = pix.child_end(sampling_level);
  for (pixel child_pix = pix.child_begin(sampling_level); child_pix
      != child_end; child_pix.next()) {
    total_area += contained_area(child_pix);
  }
  return total_area;
}

bool circle_bound::may_intersect(const pixel& pix) const {
  return get_s2cap().MayIntersect(pix.get_cell());
}

circle_bound circle_bound::get_complement() const {
  return is_empty() ? circle_bound(-axis_, 2.0) : circle_bound(-axis_, 2.0
      - height_);
}

void circle_bound::initialize_random() {
  mtrand_.seed();
  rotate_ = axis_.dot(point(0.0, 0.0, 1.0, 1.0));
  great_circle_norm_ = point(0.0, 0.0, 1.0, 1.0).cross(axis_);
}

} //end namspace s2omp
