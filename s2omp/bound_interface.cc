#include "bound_interface.h"
#include "circle_bound.h"
#include "coverer.h"
#include "pixel.h"
#include "point.h"

namespace s2omp {

bound_interface::~bound_interface() {
}

bool bound_interface::is_empty() const {
  return true;
}

long bound_interface::size() const {
  return 0;
}

void bound_interface::clear() {
}

double bound_interface::area() const {
  return 0.0;
}

bool bound_interface::contains(const point& p) const {
  return false;
}

bool bound_interface::contains(const pixel& pix) const {
  return false;
}

double bound_interface::contained_area(const pixel& pix) const {
  return 0.0;
}

bool bound_interface::may_intersect(const pixel& pix) const {
  return false;
}

point bound_interface::get_center() const {
  return point();
}

circle_bound bound_interface::get_bound() const {
  return circle_bound();
}

void bound_interface::get_covering(pixel_vector* pixels) const {
  if (!pixels->empty()) pixels->clear();
  coverer cover;

  cover.get_covering(*this, pixels);
}

void bound_interface::get_size_covering(
    long max_pixels, pixel_vector* pixels) const {
  if (!pixels->empty()) pixels->clear();
  coverer cover;

  cover.get_size_covering(max_pixels, *this, pixels);
}

void bound_interface::get_area_covering(
    double fractional_area_tolerance, pixel_vector* pixels) const {
  if (!pixels->empty()) pixels->clear();
  coverer cover;

  cover.get_area_covering(fractional_area_tolerance, *this, pixels);
}

void bound_interface::get_interior_covering(
    int max_level, pixel_vector* pixels) const {
  if (!pixels->empty()) pixels->clear();
  coverer cover;

  cover.get_interior_covering(*this, pixels);
}

// We may want to overwrite this covering over other bounds.
void bound_interface::get_simple_covering(
    int level, pixel_vector* pixels) const {
  if (!pixels->empty()) pixels->clear();
  coverer::get_simple_covering(*this, level, pixels);
}

void bound_interface::get_center_covering(
    int level, pixel_vector* pixels) const {
  if (!pixels->empty()) pixels->clear();
  coverer::get_center_covering(*this, level, pixels);
}

point bound_interface::get_random_point() {
  circle_bound cap = get_bound();

  point p = cap.create_random_point();
  while (!contains(p)) {
    p = cap.create_random_point();
  }

  return p;
}

void bound_interface::get_random_points(
    long int n_points, point_vector* points) {
  if (!points->empty()) points->clear();

  circle_bound cap = get_bound();

  for (long int i = 0; i < n_points; ++i) {
    point p = cap.create_random_point();
    while (!contains(p)) {
      p = cap.create_random_point();
    }

    points->push_back(p);
  }
}

} //end namespace s2omp
