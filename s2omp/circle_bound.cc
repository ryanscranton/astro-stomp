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
  initialized_random_ = false;
}

circle_bound::circle_bound(const point& axis, double height) {
  axis_ = axis;
  height_ = height;
  initialized_random_ = false;
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
    point p_norm = p;
    p_norm.normalize();
    double dist2 = (axis_.x() - p_norm.x()) * (axis_.x() - p_norm.x())
        + (axis_.y() - p_norm.y()) * (axis_.y() - p_norm.y()) + (axis_.z()
        - p_norm.z()) * (axis_.z() - p_norm.z());
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

  // Generate our random values following s2testing::SamplePoint
  double h = height_ * random_double();
  double theta = 2.0 * PI * random_double();
  double radius = sqrt(h * (2.0 - h)); // Radius of the circle.

  // Now we re-project the point from the axis frame to the X-Y-Z frame
  S2Point s2point = S2::FromFrame(frame_, S2Point(cos(theta) * radius, sin(
      theta) * radius, 1.0 - h));

  // And return the resulting point.
  return point::s2point_to_point(s2point.Normalize());
}

void circle_bound::get_weighted_random_points(long n_points,
    point_vector* points, const point_vector& input_points) {
  if (!points->empty())
    points->clear();

  points->reserve(n_points);
  for (long i = 0; i < n_points; ++i) {
    point p = create_random_point();

    int idx = random_int(input_points.size());
    p.set_weight(input_points[idx].weight());

    points->push_back(p);
  }
}

bool circle_bound::is_empty() const {
  return height_ < 0;
}

long circle_bound::size() const {
  return is_empty() ? 0 : 1;
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

bool circle_bound::may_intersect(const pixel& pix) const {
  return get_s2cap().MayIntersect(pix.get_cell());
}

point circle_bound::get_random_point() {
  return create_random_point();
}

void circle_bound::get_random_points(long n_points, point_vector* points) {
  if (!points->empty())
    points->clear();

  points->reserve(n_points);
  for (long i = 0; i < n_points; ++i) {
    points->push_back(create_random_point());
  }
}

circle_bound circle_bound::get_complement() const {
  return is_empty() ? circle_bound(-axis_, 2.0) : circle_bound(-axis_, 2.0
      - height_);
}

void circle_bound::initialize_random() {
  S2::GetFrame(point::point_to_s2point(axis_), &frame_);
  initialized_random_ = true;
}

} //end namspace s2omp
