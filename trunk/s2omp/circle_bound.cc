/*
 * circle_bound.cc
 *
 *  Created on: Jul 18, 2012
 *      Author: cbmorrison
 */


// #include "angular_bin-inl.h"
// #include "bound_interface.h"
#include "circle_bound.h"
#include "pixel.h"
#include "point.h"
#include "s2cell.h"

namespace s2omp {

// circle_bound::circle_bound() {
//   axis_ = point();
//   height_ = -1;
// }

circle_bound::circle_bound() {
}

circle_bound::circle_bound(const point& axis, double height) {
  axis_ = axis;
  height_ = height;
}

// TODO(cbmorrison) Commeting this out for now as we are just testing the
// core methods.
/*
circle_bound* circle_bound::from_angular_bin(const point& axis,
      const angular_bin& bin) {
  circle_bound* bound = new circle_bound(axis, 1.0 - bin.cos_theta_max());
  return bound;
}
*/

circle_bound* circle_bound::from_radius(
    const point& axis, double radius_degrees) {
  return new circle_bound(axis, 1.0 - cos(radius_degrees * DEG_TO_RAD));
}

circle_bound* circle_bound::from_height(const point& axis, double height) {
  return new circle_bound(axis, height);
}

bool circle_bound::is_empty() const {
  return height_ < 0;
}

long circle_bound::size() const {
  if (is_empty())
    return 0;
  return 1;
}

double circle_bound::area() const {
  if (!is_empty())
    return 2*PI*height_;
  return 0.0;
}

bool circle_bound::contains(const point& p) const {
  double p_height = 1.0 - axis_.dot(p);
  return height_ >= p_height;
}

bool circle_bound::contains(const pixel& pix) const {
  S2Cell cell = pix.get_cell();
  for (int k = 0; k < 4; ++k) {
    if (!contains(pixel::s2point_to_point(cell.GetVertexRaw(k)))) return false;
  }
  circle_bound* complement = get_complement();
  bool ans = !complement->may_intersect(pix);
  delete complement;
  return ans;
}

// TODO(cbmorrison) this may be a little bit slower than we want. Is there some
// way of just to make this better by returning a vector of pixels for the
// children?
double circle_bound::contained_area(const pixel& pix) const {
 double area = 0.0;
 if (contains(pix)) {
   return pix.exact_area();
 } else if (may_intersect(pix)) {
   pixel child_end = pix.child_end();
   for (pixel child_pix = pix.child_begin();
        child_pix != child_end; child_pix.next()) {
     area += contained_area(child_pix);
   }
 }
 return area;
}

// TODO(cbmorrison) the intermediate S2Cell creation may be expensive and slow.
// If needed create an s2omp class instead of wrapping cap_.may_intersect.
bool circle_bound::may_intersect(const pixel& pix) const {
  S2Cell cell = pix.get_cell();
  point_vector vertices;
  vertices.reserve(4);
  for (int k = 0; k < 4; ++k) {
    point vertex = pixel::s2point_to_point(cell.GetVertexRaw(k));
    if (contains(vertex)) {
      return true;
    } else {
      vertices.push_back(vertex);
    }
  }

  return intersects(pix, vertices);
}

point circle_bound::get_center() const {
  return axis_;
}

circle_bound circle_bound::get_bound() const {
  return circle_bound(axis_, height_);
}

point circle_bound::get_random_point() {
  if (!initialized_random_) {
    initialize_random();
  }
  // To generate a random point on the sphere we first generate the point as
  // if this cap's axis is the z axis. We then rotate the point to the true
  // position on the sphere.

  // Generate our random values. We first generate a point uniform in
  // cos(theta) from the lowest point on the cap to the highest. Next we
  // generate a point uniform in phi.
  double z = mtrand.rand(height_) + 1 - height_;
  double phi = mtrand.rand(2.0*PI);

  // To turn our angles into x,y,z coordinates we need to now the sin as well.
  double sintheta = sin(acos(z));

  // Now that we have the random point we need to rotate it to the correct
  // position on the sphere. We do this by rotating the point around the normal
  // to the great circle defined by the z axis and the cap axis. We rotate by
  // the angle the cap axis makes with the z axis.
  point p = point(sintheta*cos(phi), sintheta*sin(phi), z, 1.0);
  p.rotate_about(great_circle_norm_, rotate_);

  return p;
}

void circle_bound::get_random_points(long n_points,
    point_vector* points) {
  if (!points->empty())
    points->clear();
  points->reserve(n_points);

  for (long i = 0; i < n_points; ++i) {
    points->push_back(get_random_point());
  }
}

point circle_bound::get_weighted_random_point(
    const point_vector& points) {
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

bool circle_bound::intersects(
    const pixel& pix, const point_vector& vertices) const {
  // Much of this code is lifted from S2::S2Cap.Intersects

  // If the cap that we are considering is a hemisphere or larger, then since
  // there are no cells that would intersect this cell that would not already
  // have a vertex contained.
  if (height_ >= 1) return false;

  // Check for empty caps before checking if the axis is contained by the
  // pixel.
  if (is_empty()) return false;

  // If there are no vertices contained in the bound, but the axis is contained
  // then this bound intersects the pixel.
  if (pix.contains(axis_)) return true;

  // Since the vertices aren't contained and the axis isn't contained then
  // the only points left to test are those long the interior of an edge.

  double sin2_angle = height_ * (2 - height_);
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

circle_bound* circle_bound::get_complement() const {
  if (is_empty())
    return from_height(-axis_, 2.0);
  return from_height(-axis_, 2.0 - height_);
}

void circle_bound::initialize_random() {
  mtrand.seed();
  rotate_ = axis_.dot(point(0.0, 0.0, 1.0, 1.0));
  great_circle_norm_ = point(0.0, 0.0, 1.0, 1.0).cross(axis_);
}

} //end namspace s2omp
