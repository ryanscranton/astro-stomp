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

#include "pixel_union.h"
#include "pixel.h"
#include "point.h"


pixel_union::pixel_union() {
  min_level_ = MAX_LEVEL;
  max_level_ = 0;
  double area_ = 0.0;
  bool initialized_ = false;
}

void pixel_union::init(const pixel_vector& pixels) {
  // Most of this method is from S2CellUnion Normalize().
  if (!pixels_.empty()) {
    pixels_.clear();
  }
  pixels_.reserve(pixels.size());
  initialized_ = false;
  int min_level = MAX_LEVEL;
  int max_level = 0;
  double area = 0.0;

  // To make the process of combining and throwing out duplicate pixels we need
  // to assure that our input pixels are sorted by id.
  sort(pixels.begin(), pixels.end());

  for (pixel_iterator iter = pixels.begin(); iter != pixels.end(); ++iter) {
    // Like S2CellUnion we first check if this is a duplicate or contained
    // within a previous pixel.  If the previous pixel contains this pixel,
    // go to next position in the pixel_vector.
    if (!pixels_.empty() && pixels_.back().contains(*iter)) continue;

    // If this pixel contains any of the previous pixels, discard those.
    while (!pixels_.empty() && iter->contains(pixels_.back())) {
      pixels_.pop_back();
    }

    while (pixels_.size() >= 3) {
      // Now we begin the check to combine pixels that are children of the same
      // parent. Like in S2 we first test if the pixels are exclusive
      if ((pixels_.end()[-3].id() ^ pixels_.end()[-2].id() ^
          pixels.back().id()) != iter->id())
        break;

      // Now that we've found that the pixels are exclusive we need to perform
      // the test required o combine the pixels into one.

      // from the S2 documentation
      // "Now we do a slightly more expensive but exact test.  First, compute a
      // mask that blocks out the two bits that encode the child position of
      // "id" with respect to its parent, then check that the other three
      // children all agree with "mask."
      uint64 mask = (iter->id()).lsb() << 1;
      mask = ~(mask + (mask << 1));
      uint64 id_masked = (iter->id() & mask);
      if ((pixels_.end()[-3].id() & mask) != id_masked ||
          (pixels_.end()[-2].id() & mask) != id_masked ||
          (pixels_.end()[-1].id() & mask) != id_masked ||
          iter->is_face())
        break;

      // Now we can delete the child pixels and add the parent to the vector of
      // pixels.
      pixel_.erase(pixels_.end() -3, pixels_.end());
      *iter = iter->parent();
    }
    if (iter->level() > max_level) max_level = iter->level();
    if (iter->level() < min_level) min_level = iter->level();
    area += iter->exact_area();
    pixels_.push_back(*iter);
  }
  if (pixels_.size() < pixels.size()) {
    init_swap(&pixels_);
  }
  if (!pixels_.empty()) {
    initialized_ = true;
    max_level_ = max_level;
    min_level_ = min_level;
    area_ = area;
  }
}

void pixel_union::init_swap(pixel_vector* pixels) {
  pixels_.swap(*pixels);
  pixels->clear();
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

void pixel_union::combine(const pixel_union& u) {
  // This method is easy enough. We can simply concatenate the two pixel vectors
  // in each union and re-initialize the map. This may not be the best plan as
  // since both maps are already unions they should conform to the ordering in
  // init. Using this assumption could simplify things.

  // This method has the current problem that duplicate area may be covered by
  // pixels at different resolutions.

  pixel_vector pixels;
  pixels.reserve(size() + u.size());
  pixel_iterator iter = begin();
  pixel_iterator u_iter = u.begin();

  while (iter != end() || u_iter != u.end()) {
    if (u_iter == u.end || iter->id() < u_iter->id()) {
      pixels.push_back(*iter);
      ++iter;
      continue;
    } else if (iter != end || iter->id() > u_iter->id()) {
      pixels.push_back(*u_iter);
      ++u_iter;
      continue;
    } else if (iter != end && u_iter != u.end && iter->id() == u_iter->id()) {
      pixels.push_back(*iter);
      ++iter;
      ++u_iter;
      continue;
    }
  }
  init(pixels);
}

void pixel_union::intersect(const pixel_union& u) {
  // loop through of input union testing first for may_interect with each this
  // union. The next test is contains and if it is then, we keep the pixel. If
  // the test is only intersects then we break the pixel up and perform the
  // same test on it's children until contains is true or intersects is
  // false.

  pixel_vector pixels;
  pixel_iterator iter = begin();
  pixel_iterator u_iter = u.being();

  while (iter != end() || u_iter != u.end()) {
    if (u_iter == u.end() || iter->id() < u_iter->id()) {
      if (u.contains(*iter)) pixels.push_back(*iter);
      ++iter;
      continue;
    } else if (iter == end() || iter->id() > u_iter->id()) {
      if (contains(*u_iter)) pixels.push_back(*u_iter);
      ++u_iter;
      continue;
    } else if (iter != end() && u_iter != u.end && iter->id() == u_iter->id()) {
      pixels.push_back(*iter);
      ++iter;
      ++u_iter;
      continue;
    }
  }
  init(pixels);
}


void pixel_union::exclude(const pixel_union& u) {
  // same as above but we want to flip the tests so that if a pixel in our
  // current union is within the input union we want to drop that pixel.

  pixel_vector* pixels;
  for (pixel_iterator iter = begin(); iter != end(); ++iter) {
    pixel_exclusion(*iter, u, pixels);
  }
  init(*pixels);
}

void pixel_union::init_from_combination(const pixel_union& a,
    const pixel_union& b) {

  pixel_vector pixels;
  pixels.reserve(a.size() + b.size());
  pixel_iterator a_iter = a.begin();
  pixel_iterator b_iter = b.begin();

  while (a_iter != a.end() || b_iter != b.end()) {
    if (b_iter == b.end || a_iter->id() < b_iter->id()) {
      pixels.push_back(*a_iter);
      ++iter;
      continue;
    } else if (a_iter == a.end() || a_iter->id() > b_iter->id()) {
      pixels.push_back(*b_iter);
      ++b_iter;
      continue;
    } else if (a_iter != end && b_iter != b.end() &&
        a_iter->id() == b_iter->id()) {
      pixels.push_back(*a_iter);
      ++a_iter;
      ++b_iter;
      continue;
    }
  }
  init(pixels);
}

void pixel_union::init_from_intersection(const pixel_union& a,
    const pixel_union& b) {

  pixel_vector pixels;
  pixel_iterator a_iter = a.begin();
  pixel_iterator b_iter = u.begin();

  while (a_iter != a.end() || b_iter != b.end()) {
    if (b_iter == b.end() || a_iter->id() < b_iter->id()) {
      if (b.contains(*a_iter)) pixels.push_back(*a_iter);
      ++a_iter;
      continue;
    } else if (a_iter == end() || a_iter->id() > b_iter->id()) {
      if (a.contains(*b_iter)) pixels.push_back(*b_iter);
      ++b_iter;
      continue;
    } else if (a_iter != end() && b_iter != u.end &&
        a_iter->id() == b_iter->id()) {
      pixels.push_back(*a_iter);
      ++a_iter;
      ++b_iter;
      continue;
    }
  }
  init(pixels);
}

void pixel_union::init_from_exclusion(const pixel_union& a,
    const pixel_union& b) {
  pixel_vector* pixels;
  for (pixel_iterator iter = begin(); iter != end(); ++iter) {
    a.pixel_exclusion(*iter, &b, pixels);
  }
  init(*pixels);
}

bool pixel_union::contains(point& p) {
  return contains(p.to_pixel(MAX_LEVEL));
}

bool pixel_union::contains(pixel& pix) {

  // again this is already done pretty elegantly in S2 so we'll just adapt it
  // for stomp.
  pixel_iterator iter = lower_bound(begin(), end(), pix);
  if (iter != end() && iter->range_min() <= pix) return true;
  return iter != begin() && (--iter)->range_max() >= pix;
}

bool pixel_union::intersects(pixel& pix) {

  pixel_iterator iter = lower_bound(begin(), end(), pix);
  if (iter != end() && iter->range_min() <= pix.range_max()) return true;
  return iter != begin() && (--iter)->range_max() >= pix.range_min();
}

void pixel_union::initialzie_bound() {
  if (pixels_.empty()) return cirlce_bound();
  point center(0, 0, 0, 1);
  for (pixel_iterator iter = begin(); iter != end(); ++iter) {
    center += pixel.exact_area()*iter->to_point();
  }
  if (center == point(0, 0, 0, 1)) {
    center = point(1, 0, 0);
  } else {
    center = center.norm();
  }

  circle_bound bound = cirlce_bound.from_axis_height(center, 0);
  for (pixel_iterator iter = begin(); iter != end(); ++iter) {
    bound.add_cap(iter->get_bound());
  }
  bound_ = bound;
}

point pixel_union::generate_random_point() const {
  if (!initialized_bound_) initialze_bound();

  point p = bound_.generate_random_point();
  while (!contains(p)) {
    p = bound_.generate_random_point();
  }
  return p;
}

void pixel_union::generate_random_points(long int n_points,
    points_vector* points) const {
  if (!points->empty()) points->clear();

  for (long int i = 0; i < n_points; ++i) {
    points->push_back(generate_random_point());
  }
}

void pixel_union::generate_weighted_random_point(const point_vector& points) {
  if (!initialized_bound_) initialze_bound();

  point p = bound_.generate_weighted_random_point(points);
  while (!contains(p)) {
    p = bound_.generate_weighted_random_point(points);
  }
  return p;
}

void pixel_union::generate_weighted_random_points(long n_points,
    const point_vector& points, point_vector* points) {
  if (!points->empty()) points->clear();

  for (long int i = 0; i < n_points; ++i) {
    points->push_back(generate_weighted_random_point());
  }
}

virtual void covering(pixel_vector* pixels) const {
  // this is the default mode for a pixel covering. For this as for other
  // coverings we use only 8 pixels at max to cover our union.
  covering(8, pixels);
}

virtual void covering(int max_pixels, pixel_vector* pixels) const {
  // For this class we want to keep as few pixels as possible (defined by
  // max_pixels) and retain a close approximation of the area contained by the
  // union.
  if (!pixels->empty()) pixels->clear();
  double average_area = area_/(1.0*max_pixels);
  int level = MAX_LEVEL;
  double pixel_average = -1.0;

  while (pixel.average_area(level) < average_area) {
    if (level == 0) break;
    level--;
  }
  if (level < 30) level++;

  while (pixels->empty() || pixels->size() > max_pixels) {
    pixels->clear();
    simple_covering(level, pixels);
    level--;
  }
}

void pixel_union::simple_covering(int level, pixel_vector* pixels) const {
  // For a simple covering of a pixel_union we don't need anything fancy like
  // FloodFill from S2 (at least for this method). If we want a simple covering
  // of a pixel_union then we just need to loop through the pixels, make parents
  // out of children
  if (!pixels->empty()) pixels->clear();

  for (pixel_iterator iter = begin(); iter != end(); ++iter) {
    if (iter->level() < level) {
      for (pixel_iterator child_iter = child_iter->child_begin(level);
          child_iter != child_end; ++child_iter) {
        if (!pixels->empty() && pixels->back().contains(*child_iter)) continue;
        pixels->push_back(*child_iter);
      }
    } else if (iter->level() == level) {
      if (!pixels->empty() && pixels->back().contains(*iter)) continue;
      pixels->push_back(*iter);
    } else if (iter->level() > level) {
      pixel par = iter->parent(level);
      if (!pixels->empty() && pixels->back().contains(par)) continue;
      pixels->push_back(par);
    }
  }
}

void pixel_union::pixel_intersection(const pixel& pix, const pixel_union& u,
    pixel_vector* pixels) {
  if (u->intersects(pix)) {
    if (u->contains(pix)) {
      pixels->push_back(pix);
    } else if (!pix.is_leaf()) {
      for (pixel_iterator child_iter = pix.child_begin(); child_iter
          != pix.child_end(); ++child_iter) {
        pixel_intersection(*child_iter, u, pixels);
      }
    }
  }
}

void pixel_union::pixel_exclusion(const pixel& pix, const pixel_union& u,
    pixel_vector* pixels) {
  if (!u->contains(pix)) {
    if (!u->intersects(pix)) {
      pixels->push_back(pix);
    } else if (!pix.is_leaf()) {
      for (pixel_iterator child_iter = pix.child_begin(); child_iter
          != pix.child_end(); ++child_iter) {
        pixel_exclusion(*child_iter, u, pixels);
      }
    }
  }
}

virtual bool may_intersect(const pixel& pix) const {
  bool intersected = false;
  for (pixel_iterator iter = begin(); iter != end(); ++iter)
    if (iter->may_intersect())
      intersected = true;
  return intersected;
}
