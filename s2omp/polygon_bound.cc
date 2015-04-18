#include <s2/s2cell.h>
#include <s2/s2loop.h>
#include <s2/s2polygon.h>

#include "polygon_bound.h"

#include "angular_bin-inl.h"
#include "circle_bound.h"
#include "pixel.h"
#include "point.h"

namespace s2omp {

polygon_bound::polygon_bound() {
  polygon_ = new S2Polygon();
}

polygon_bound::polygon_bound(const point_vector& points) {
  if (points.empty()) {
    polygon_ = new S2Polygon();
    return;
  }

  std::vector<S2Loop*> loops;
  loops.push_back(polygon_bound::point_vector_to_s2loop(points));
  polygon_ = new S2Polygon(&loops);
}

polygon_bound::~polygon_bound() {
  delete polygon_;
}

polygon_bound* polygon_bound::from_points(const point_vector& points) {
  return new polygon_bound(points);
}

bool polygon_bound::add_loop(const point_vector& points) {
  std::vector<S2Loop*> loops;
  // Copy the polygon loops to our vector and reset the polygon.
  polygon_->Release(&loops);

  // Add our new loop to the list and re-initialize.
  loops.push_back(polygon_bound::point_vector_to_s2loop(points));
  polygon_->Init(&loops);

  // Return true if we have a valid polygon.
  return polygon_->IsValid();
}

bool polygon_bound::intersects(const polygon_bound& bound) const {
  return polygon_->Intersects(bound.get_s2polygon());
}

S2Loop* polygon_bound::point_vector_to_s2loop(const point_vector& points) {
  std::vector<S2Point> s2points;
  s2points.reserve(points.size());
  for (int i = 0; i < points.size(); i++) {
    s2points.push_back(point::point_to_s2point(points[i]).Normalize());
  }

  return new S2Loop(s2points);
}

bool polygon_bound::is_empty() const {
  return polygon_->num_vertices() == 0;
}

long polygon_bound::size() const {
  return polygon_->num_vertices();
}

void polygon_bound::clear() {
  delete polygon_;
  polygon_ = new S2Polygon();
}

double polygon_bound::area() const {
  return polygon_->GetArea() * STRAD_TO_DEG2;
}

bool polygon_bound::contains(const point& p) const {
  return polygon_->Contains(point::point_to_s2point(p));
}

bool polygon_bound::contains(const pixel& pix) const {
  return polygon_->Contains(pix.get_cell());
}

double polygon_bound::contained_area(const pixel& pix) const {
  S2Cell cell = pix.get_cell();
  if (polygon_->Contains(cell)) {
    return pix.exact_area();
  } else if (!polygon_->MayIntersect(cell)) {
    return 0.0;
  }

  S2Polygon* cell_polygon = new S2Polygon(cell);
  S2Polygon* intersection = new S2Polygon();
  intersection->InitToIntersection(polygon_, cell_polygon);

  double area = intersection->GetArea() * STRAD_TO_DEG2;
  delete cell_polygon;
  delete intersection;

  return area;
}

bool polygon_bound::may_intersect(const pixel& pix) const {
  return polygon_->MayIntersect(pix.get_cell());
}

point polygon_bound::get_center() const {
  return point::s2point_to_point(polygon_->GetCentroid());
}

circle_bound polygon_bound::get_bound() const {
  S2Cap cap = polygon_->GetCapBound();
  return circle_bound(point::s2point_to_point(cap.axis()), cap.height());
}

} // end namespace s2omp
