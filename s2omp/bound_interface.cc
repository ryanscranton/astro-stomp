/*
 * bound_interface.cc
 *
 *  Created on: Aug 17, 2012
 *      Author: morrison
 */


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

void bound_interface::get_covering(const uint32_t max_pixels,
    pixel_vector* pixels) const {
  if (!pixels->empty()) pixels->clear();
  coverer cover;

  cover.get_covering(max_pixels, *this, pixels);
}

void bound_interface::get_covering(double fractional_area_tolerance,
    pixel_vector* pixels) const {
  if (!pixels->empty()) pixels->clear();
  coverer cover;

  cover.get_covering(fractional_area_tolerance, *this, pixels);
}

void bound_interface::get_interior_covering(int max_level,
    pixel_vector* pixels) const {
  if (!pixels->empty()) pixels->clear();
  coverer cover;

  cover.get_interior_covering(*this, pixels);
}

// We may want to overwrite this covering over other bounds.
void bound_interface::get_simple_covering(
    int level, pixel_vector* pixels) const {
  if (!pixels->empty()) pixels->clear();
  coverer cover;

  cover.get_simple_covering(*this, level, pixels);
}

} //end namespace s2omp
