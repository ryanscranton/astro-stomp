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

void bound_interface::get_covering(pixel_vector* pixels) {
  if (!pixels->empty()) pixels->clear();
  coverer cover = coverer();

  cover.get_covering(*this, pixels);
}

void bound_interface::get_covering(const uint32_t max_pixels,
    pixel_vector* pixels) {
  if (!pixels->empty()) pixels->clear();
  coverer cover = coverer();

  cover.get_covering(max_pixels, *this, pixels);
}

void bound_interface::get_covering(double fractional_area_tolerance,
    pixel_vector* pixels) {
  if (!pixels->empty()) pixels->clear();
  coverer cover = coverer();

  cover.get_covering(fractional_area_tolerance, *this, pixels);
}

void bound_interface::get_interior_covering(int max_level,
    pixel_vector* pixels) {
  if (!pixels->empty()) pixels->clear();
  coverer cover = coverer(0, max_level);

  cover.get_interior_covering(*this, pixels);
}

// We may want to overwrite this covering over other bounds.
void bound_interface::get_simple_covering(int level,
    pixel_vector* pixels) {
  if (!pixels->empty()) pixels->clear();
  coverer cover = coverer();

  cover.get_simple_covering(*this, level, pixels);
}

} //end namespace s2omp
