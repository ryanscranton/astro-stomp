/*
 * circle_bound.cc
 *
 *  Created on: Jul 18, 2012
 *      Author: cbmorrison
 */




#include "circle_bound.h"

circle_bound::circle_bound() {
  radius_ = -1.0;
  height_ = 100;
  initialized_random_ = false;
}

circle_bound* circle_bound::from_radius(point& axis, double radius_degrees) {
}

bool circle_bound::contains(const point& p) {
  return axis_.dot(p) >= height_;
}

bool circle_bound::contains(const pixel& pix) {
  // maybe check all the vertices an that is sufficient?
  // or return the cap associated with the pixel and check contains?
}

double circle_bound::contained_area(const pixel& px) {
 double area = 0.0;
 if (contains(pix)) {
   area += pix.exact_area();
 } else if (may_intersect(pix)) {
   for (pixel_iterator iter = pix.child_begin();
       iter != pix.child_end(); ++iter) {
     area += contained_area(*iter);
   }
 }
 return area;
}

point circle_bound::generate_random_point() {
  if (!initialized_random_) {
    initialize_random();
  }
  // To generate a random point on the sphere we first generate the point as
  // if this cap's axis is the z axis. We then rotate the point to the true
  // position on the sphere.

  // Generate our random values. We first generate a point uniform in
  // cos(theta) from the lowest point on the cap to the highest. Next we
  // generate a point uniform in phi.
  double z = mtrand.rand(1.0 - height_) + height_;
  double phi = mtrand.rand(2.0*PI);

  // To turn our angles into x,y,z coordinates we need to now the sin as well.
  double sintheta = sin(acos(z));

  // Now that we have the random point we need to rotate it to the correct
  // position on the sphere. We do this by rotating the point around the normal
  // to the great circle defined by the z axis and the cap axis. We rotate by
  // the angle the cap axis makes with the z axis.
  point p = point(sintheta*cos(phi), sintheta*sin(phi), z);
  p.rotate_about(great_norm_, rotate_);

  return p;
}

void cirlce_bound::generate_random_points(long n_points, point_vector* points) {
  if (!points->empty()) points->clear();
  for (long i = 0; i < n_points; ++i) {
    points->push_back(generate_random_point());
  }
}
