// Copyright 2012  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)
//         cmmorrison@gmail.com (Chris Morrison)
//
// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains the Map class.  Maps are intended to describe an
// arbitrary area on the sky as precisely as possible, given the limits of
// pixel resolution, physical memory, etc.  Internally, this information is
// encoded using pixels at a range of resolutions, depending on how the
// boundaries of the area in question and the pixelization scheme interact.
// However, the goal of the class is to abstract away those details, allowing
// the user to treat Maps as a pure representative of spherical geometry.

#include "core.h"
#include "pixel_union.h"
#include "pixel.h"
#include "point.h"

#include <s2/s2cellunion.h>

namespace s2omp {

pixel_union::pixel_union() {
  min_level_ = MAX_LEVEL;
  max_level_ = 0;
  area_ = 0.0;
  initialized_ = false;
}

pixel_union::~pixel_union() {
  pixels_.clear();
  bound_.clear();
  initialized_ = false;
  min_level_ = MAX_LEVEL;
  max_level_ = 0;
  area_ = 0.0;
}

void pixel_union::init(pixel_vector& pixels) {
  clear();
  if (pixels.empty()) {
    return;
  }

  std::vector<S2CellId> input_cells;
  std::vector<S2CellId> output_cells;
  input_cells.reserve(pixels.size());
  output_cells.reserve(pixels.size());

  // Copy the S2CellIds into a vector and initialize an S2CellUnion from them.
  for (pixel_iterator iter = pixels.begin(); iter != pixels.end(); ++iter) {
    input_cells.push_back(iter->get_cellid());
  }
  pixels.clear();

  S2CellUnion s2cellunion;
  s2cellunion.Init(input_cells);
  input_cells.clear();

  // Now extract the normalized S2CellIds and convert them to pixels.
  s2cellunion.Detach(&output_cells);

  pixels_.reserve(output_cells.size());
  initialized_ = false;
  int min_level = MAX_LEVEL;
  int max_level = 0;
  double area = 0.0;
  for (int k = 0; k < output_cells.size(); k++) {
    pixel pix(output_cells[k]);
    if (pix.level() > max_level)
      max_level = pix.level();
    if (pix.level() < min_level)
      min_level = pix.level();
    area += pix.exact_area();
    pixels_.push_back(pix);
  }

  if (!pixels_.empty()) {
    initialized_ = true;
    max_level_ = max_level;
    min_level_ = min_level;
    area_ = area;
    range_min_ = pixels_.front().range_min();
    range_max_ = pixels_.back().range_max();
  }
}

void pixel_union::soften(int max_level) {
  // Since init has already given us an ordered set of pixels all we need to do
  // is loop through the pixels until we hit the resolution we are interested in
  // softening to. Starting from the end of the vector and looping backwards
  // sounds like the best bet here.

  pixel_vector pixels;
  for (pixel_iterator iter = pixels_.begin(); iter != pixels_.end(); ++iter) {
    if (iter->level() <= max_level) {
      pixels.push_back(*iter);
    } else {
      if (!pixels.empty() && pixels.back().contains(iter->parent(max_level))) {
        continue;
      }
      pixels.push_back(iter->parent(max_level));
    }
  }
  init(pixels);
}

void pixel_union::combine_with(const pixel_union& pix_union) {
  init_from_combination(*this, pix_union);
}

void pixel_union::intersect_with(const pixel_union& pix_union) {
  init_from_intersection(*this, pix_union);
}

void pixel_union::exclude_from(const pixel_union& pix_union) {
  init_from_exclusion(*this, pix_union);
}

void pixel_union::init_from_combination(const pixel_union& a,
    const pixel_union& b) {
  // We could do some work here to merge sort the two vector sets, but since
  // the initialization will sort them anyway, that effort's likely wasted.
  // Instead, just concatenate the two pixel_vectors and initialize.
  pixel_vector pixels;
  pixels.reserve(a.size() + b.size());
  pixels.insert(pixels.end(), a.begin(), a.end());
  pixels.insert(pixels.end(), b.begin(), b.end());
  init(pixels);
}

void pixel_union::init_from_intersection(const pixel_union& a,
    const pixel_union& b) {
  // If there's definitely no intersection between the two unions, then stop.
  if (!a.may_intersect(b)) {
    clear();
    return;
  }

  // If there's possibly some intersection, then we try to find it.
  pixel_vector pixels;
  if (a.size() < b.size()) {
    for (pixel_iterator iter = b.begin(); iter != b.end(); ++iter) {
      pixel_vector temp_pixels;
      a.pixel_intersection(*iter, &temp_pixels);
      pixels.insert(pixels.end(), temp_pixels.begin(), temp_pixels.end());
    }
  } else {
    for (pixel_iterator iter = a.begin(); iter != a.end(); ++iter) {
      pixel_vector temp_pixels;
      b.pixel_intersection(*iter, &temp_pixels);
      pixels.insert(pixels.end(), temp_pixels.begin(), temp_pixels.end());
    }
  }

  init(pixels);
}

void pixel_union::init_from_exclusion(const pixel_union& a,
    const pixel_union& b) {
  pixel_vector pixels;
  // Unlike combination or intersection, order is important here.  We want the
  // part of A that is outside of B, not vice versa.
  for (pixel_iterator iter = a.begin(); iter != a.end(); ++iter) {
    pixel_vector temp_pixels;
    b.pixel_exclusion(*iter, &temp_pixels);
    pixels.insert(pixels.end(), temp_pixels.begin(), temp_pixels.end());
  }

  init(pixels);
}

bool pixel_union::intersects(const pixel& pix) const {
  pixel_iterator iter = lower_bound(begin(), end(), pix);
  if (iter != end() && iter->range_min() <= pix.range_max())
    return true;
  return iter != begin() && (--iter)->range_max() >= pix.range_min();
}

bool pixel_union::may_intersect(const pixel_union& pix_union) const {
  // If the range_min and range_max for the two union overlap, then they may
  // intersect.  This is a faster check than verifying intersection.
  return (range_min_ <= pix_union.range_min() && range_max_
      >= pix_union.range_min()) || (range_min_ <= pix_union.range_max()
      && range_max_ >= pix_union.range_max());
}

bool pixel_union::intersects(const pixel_union& pix_union) const {
  if (!may_intersect(pix_union)) {
    return false;
  }

  for (pixel_iterator iter = pixels_.begin(); iter != pixels_.end(); ++iter) {
    if (pix_union.intersects(*iter)) {
      return true;
    }
  }

  return false;
}

void pixel_union::pixel_intersection(const pixel_union& pix_union,
    const pixel& pix, pixel_vector* pixels) {
  if (!pixels->empty()) {
    pixels->clear();
  }

  // If the pixel doesn't intersect our union, we're done.
  if (!pix_union.intersects(pix)) {
    return;
  }

  // If the pixel is contained by the union, then add it to the output pixel
  // and we're done.
  if (pix_union.contains(pix)) {
    pixels->push_back(pix);
    return;
  }

  // If the pixel isn't fully contained, but does intersect, then we need to
  // check the children, provided that this isn't a leaf cell.
  if (pix.is_leaf()) {
    return;
  }

  for (pixel c = pix.child_begin(); c != pix.child_end(); c = c.next()) {
    // Need a temporary pixel_vector since the input pixel_vector is cleared.
    pixel_vector temp_pixels;
    pixel_intersection(pix_union, c, &temp_pixels);
    pixels->insert(pixels->end(), temp_pixels.begin(), temp_pixels.end());
  }
}

void pixel_union::pixel_intersection(const pixel& pix, pixel_vector* pixels) const {
  pixel_intersection(*this, pix, pixels);
}

void pixel_union::pixel_exclusion(const pixel_union& pix_union,
    const pixel& pix, pixel_vector* pixels) {
  if (!pixels->empty()) {
    pixels->clear();
  }

  // If the input pix is fully contained by our union, then we're done.
  if (pix_union.contains(pix)) {
    return;
  }

  // If there's no intersection between pixel and union, then we add the pixel
  // to the output pixel_vector and we're finished.
  if (!pix_union.intersects(pix)) {
    pixels->push_back(pix);
    return;
  }

  // If there is some intersection but not full containment, then we need to
  // check the child pixels, provided that this isn't a leaf cell.
  if (pix.is_leaf()) {
    return;
  }

  for (pixel c = pix.child_begin(); c != pix.child_end(); c = c.next()) {
    // Need a temporary vector here since the input pixel_vector is cleared.
    pixel_vector temp_pixels;
    pixel_exclusion(pix_union, c, &temp_pixels);
    pixels->insert(pixels->end(), temp_pixels.begin(), temp_pixels.end());
  }
}

void pixel_union::pixel_exclusion(const pixel& pix, pixel_vector* pixels) const {
  pixel_exclusion(*this, pix, pixels);
}

void pixel_union::initialize_bound() {
  if (pixels_.empty()) {
    bound_ = circle_bound();
    return;
  }

  bound_ = get_bound();
  initialized_bound_ = true;
}

void pixel_union::clear() {
  pixels_.clear();
  bound_.clear();
  min_level_ = MAX_LEVEL;
  max_level_ = 0;
  area_ = 0.0;
  initialized_ = false;
  initialized_bound_ = false;
}

bool pixel_union::contains(const point& p) const {
  return contains(p.to_pixel());
}

bool pixel_union::contains(const pixel& pix) const {
  // again this is already done pretty elegantly in S2 so we'll just adapt it
  // for stomp.
  pixel_iterator iter = lower_bound(pixels_.begin(), pixels_.end(), pix);
  if (iter != pixels_.end() && iter->range_min() <= pix)
    return true;
  return iter != pixels_.begin() && (--iter)->range_max() >= pix;
}

double pixel_union::contained_area(const pixel& pix) const {
  if (!intersects(pix)) {
    return 0.0;
  }

  pixel_vector pixels;
  pixel_intersection(pix, &pixels);
  double area = 0.0;
  for (pixel_iterator iter = pixels.begin(); iter != pixels.end(); ++iter) {
    area += iter->exact_area();
  }

  return area;
}

bool pixel_union::may_intersect(const pixel& pix) const {
  return intersects(pix);
}

point pixel_union::get_center() const {
  if (initialized_bound_) {
    return bound_.axis();
  }

  double x = 0.0, y = 0.0, z = 0.0;
  for (pixel_iterator iter = begin(); iter != end(); ++iter) {
    double area = iter->exact_area();
    point center = iter->center_point();
    x += area * center.unit_sphere_x();
    y += area * center.unit_sphere_y();
    z += area * center.unit_sphere_z();
  }

  return point(x, y, z, 1.0);
}

circle_bound pixel_union::get_bound() const {
  if (initialized_bound_) {
    return bound_;
  }

  circle_bound bound = circle_bound(get_center(), 0.0);
  for (pixel_iterator iter = begin(); iter != end(); ++iter) {
    bound.add_circle_bound(iter->get_bound());
  }

  return bound;
}

void pixel_union::get_covering(pixel_vector* pixels) const {
  // this is the default mode for a pixel covering. For this as for other
  // coverings we use only 8 pixels at max to cover our union.
  get_covering(DEFAULT_COVERING_PIXELS, pixels);
}

void pixel_union::get_covering(long max_pixels, pixel_vector* pixels) const {
  // For this class we want to keep as few pixels as possible (defined by
  // max_pixels) and retain a close approximation of the area contained by the
  // union.
  if (!pixels->empty())
    pixels->clear();
  double average_area = area_ / (1.0 * max_pixels);
  int level = MAX_LEVEL;

  while (pixel::average_area(level) < average_area) {
    if (level == 0)
      break;
    level--;
  }
  if (level < MAX_LEVEL)
    level++;

  while (pixels->empty() || pixels->size() > max_pixels) {
    pixels->clear();
    get_simple_covering(level, pixels);
    level--;
  }
}

void pixel_union::get_covering(double fractional_area_tolerance,
    pixel_vector* pixels) const {
  if (!pixels->empty())
    pixels->clear();
  pixels->reserve(pixels_.size());

  for (pixel_iterator iter = pixels_.begin(); iter != pixels_.end(); ++iter) {
    pixels->push_back(*iter);
  }
}

void pixel_union::get_interior_covering(int max_level, pixel_vector* pixels) const {
  if (!pixels->empty())
    pixels->clear();
  pixels->reserve(pixels_.size());

  for (pixel_iterator iter = pixels_.begin(); iter != pixels_.end(); ++iter) {
    if (iter->level() <= max_level) {
      pixels->push_back(*iter);
    }
  }
}

void pixel_union::get_simple_covering(int level, pixel_vector* pixels) const {
  // For a simple covering of a pixel_union we don't need anything fancy like
  // FloodFill from S2 (at least for this method). If we want a simple covering
  // of a pixel_union then we just need to loop through the pixels, make parents
  // out of children
  if (!pixels->empty())
    pixels->clear();

  for (pixel_iterator iter = pixels_.begin(); iter != pixels_.end(); ++iter) {
    if (iter->level() < level) {
      for (pixel c = iter->child_begin(); c != iter->child_end(); c = c.next()) {
        if (!pixels->empty() && pixels->back().contains(c))
          continue;
        pixels->push_back(c);
      }
    } else {
      pixel pix = iter->parent(level);
      if (!pixels->empty() && pixels->back().contains(pix))
        continue;
      pixels->push_back(pix);
    }
  }
}

void pixel_union::get_center_covering(int level, pixel_vector* pixels) const {
  get_simple_covering(level, pixels);
}

} // end namespace s2omp