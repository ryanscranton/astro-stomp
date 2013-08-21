/*
 * point.cc
 *
 *  Created on: Aug 16, 2012
 *      Author: morrison
 */

#include "pixel.h"
#include "point.h"

namespace s2omp {

point::point() {
  point_ = S2Point(0, 0, 0);
  weight_ = 0;
}

point::point(double x, double y, double z, double weight) {
  weight_ = weight;
  double norm = sqrt(x*x + y*y + z*z);
  point_ = S2Point(x/norm, y/norm, z/norm);
}

point::point(S2Point point, double weight) {
  point_ = point;
  weight_ = weight;
}

point* point::copy_point(const point& p) {
  return new point(p.s2point(), p.weight());
}

point point::copy_point(point* p) {
  return point(p->s2point(), p->weight());
}

double point::dot(const point& a, const point& b) {
  return a.unit_sphere_x() * b.unit_sphere_x() + a.unit_sphere_y()
      * b.unit_sphere_y() + a.unit_sphere_z() * b.unit_sphere_z();
}

double point::angular_distance(const point& p) const {
  return acos(dot(p));
}

double point::angular_distance(const point& a, const point& b) {
  return acos(dot(a, b));
}

// TODO(cbmorrison) Not sure what the rules should be for returning the weight
// on the cross product of the points. Currently I return the weight of the
// first point.
point point::cross(const point& p) const {
  double x = point_.y()*p.unit_sphere_z() - point_.z()*p.unit_sphere_y();
  double y = point_.z()*p.unit_sphere_x() - point_.x()*p.unit_sphere_z();
  double z = point_.x()*p.unit_sphere_y() - point_.y()*p.unit_sphere_x();

  return point(x, y, z, weight_);
}

point point::cross(const point& a, const point& b) {
  double x = a.unit_sphere_y()*b.unit_sphere_z() -
      a.unit_sphere_z()*b.unit_sphere_y();
  double y = -(a.unit_sphere_x()*b.unit_sphere_z() -
               a.unit_sphere_z()*b.unit_sphere_x());
  double z = a.unit_sphere_x()*b.unit_sphere_y() -
      a.unit_sphere_y()*b.unit_sphere_x();

  return point(x, y, z, a.weight());
}

// TODO(cbmorrison) Currently great_circle is just a wrapper for cross that
// sets the weight of the returned point to zero.
point point::great_circle(const point& p) {
  point great = cross(p);
  great.set_weight(0);
  return cross(great);
}

point point::great_circle(const point& a, const point& b) {
  point great = cross(a, b);
  great.set_weight(0);
  return great;
}

pixel point::to_pixel() const {
  return pixel(id());
}

pixel point::to_pixel(int level) const {
  return pixel(id(level));
}

void point::rotate_about(const point& axis, double rotation_angle_degrees) {
	// Time to find that python code I wrote so I can make this easy.
}




} // end namespace s2omp
