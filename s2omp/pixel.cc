// Copyright 2012  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)
//         cbmorrison@gmail.com (Chris Morrison)
//
// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains the central class for the library: the hierarchical
// pixelization that makes all of the rest of the spatial classes work.  This
// class defines all of the core Pixel operations: converting angular position
// on the sphere into pixel index, increasing and decreasing pixel resolution,
// finding neighboring pixels and so on.

#include "circle_bound.h"
#include "pixel.h"
#include "point.h"

namespace s2omp {

pixel::pixel() {
  id_ = S2CellId(0);
}

pixel pixel::from_point(const point& p) {
  return pixel(S2CellId::FromPoint(p.s2point()));
}

pixel pixel::from_point(const point& p, int level) {
  return from_point(p).parent(level);
}

pixel pixel::parent() const {
  return pixel(id_.parent());
}

pixel pixel::parent(int level) const {
  return pixel(id_.parent(level));
}

void pixel::children(pixel_vector* child_pixels) const {
  child_pixels->clear();
  if (!is_leaf()) {
    child_pixels->reserve(4);
    for (int k = 0; k < 4; k++) {
      child_pixels->push_back(pixel(id_.child(k)));
    }
  }
}

void pixel::children(int level, pixel_vector* child_pixels) const {
  child_pixels->clear();
  if (!is_leaf() && level <= MAX_LEVEL) {
    for (pixel c = child_begin(level); c != child_end(level); c = c.next()) {
      child_pixels->push_back(c);
    }
  }
}

pixel pixel::child_begin() const {
  return pixel(id_.child_begin());
}

pixel pixel::child_begin(int level) const {
  return pixel(id_.child_begin(level));
}

pixel pixel::child_end() const {
  return pixel(id_.child_end());
}

pixel pixel::child_end(int level) const {
  return pixel(id_.child_end(level));
}

pixel pixel::next() const {
  return pixel(id_.next());
}

pixel pixel::prev() const {
  return pixel(id_.prev());
}

pixel pixel::next_wrap() const {
  return pixel(id_.next_wrap());
}

pixel pixel::prev_wrap() const {
  return pixel(id_.prev_wrap());
}

bool pixel::is_cohort(const pixel& pix) const {
  return are_cohorts(*this, pix);
}

bool pixel::are_cohorts(const pixel& pix_a, const pixel& pix_b) {
  // Borrowing this implementation from S2::S2CellIdUnion::Normalize

  uint64 mask = pix_b.lsb() << 1;
  mask = ~(mask + (mask << 1));
  uint64 id_masked = (pix_b.id() & mask);
  return (pix_a.id() & mask) == id_masked;
}

bool pixel::are_cohorts(const pixel& pix_a, const pixel& pix_b,
    const pixel& pix_c, const pixel& pix_d) {
  // Borrowing this implementation from S2::S2CellIdUnion::Normalize

  // First do some fast XOR math to see if the current pixel is possibly the
  // 4th child pixel.  This establishes that the pixels are consecutive
  // bit values, which is a necessary but not sufficient condition.
  if ((pix_a.id() ^ pix_b.id() ^ pix_c.id()) != pix_d.id()) {
    return false;
  }

  // Ok, if we passed that filter, then we need a more precise but
  // expensive test.  For this, we need to create a mask that strips out
  // the bits other than the ones that correspond to our current cohort.
  // If all of our potential cohorts agree with that mask, then we have
  // a full set of children.
  uint64 mask = pix_d.lsb() << 1;
  mask = ~(mask + (mask << 1));
  uint64 id_masked = (pix_d.id() & mask);
  return (pix_a.id() & mask) == id_masked && (pix_b.id() & mask) == id_masked
      && (pix_c.id() & mask) == id_masked;
}

S2Cell pixel::get_cell() const {
  return S2Cell(id_);
}

double pixel::average_area() const {
  return get_cell().AverageArea() * STRAD_TO_DEG2;
}

double pixel::exact_area() const {
  return get_cell().ApproxArea() * STRAD_TO_DEG2;
}

bool pixel::contains(const point& p) const {
  return id_.contains(S2CellId::FromPoint(p.s2point()));
}

bool pixel::contains(const pixel& pix) const {
  return id_.contains(pix.get_cellid());
}

bool pixel::intersects(const pixel& pix) const {
  return id_.intersects(pix.get_cellid());
}

pixel pixel::range_min() const {
  return pixel(id_.range_min());
}

pixel pixel::range_max() const {
  return pixel(id_.range_max());
}

point pixel::center_point() const {
  return point::s2point_to_point(id_.ToPoint());
}

point pixel::vertex(int k) const {
  return point::s2point_to_point(get_cell().GetVertexRaw(k));
}

point pixel::edge(int k) const {
  return point::s2point_to_point(get_cell().GetEdgeRaw(k));
}

bool pixel::edge_distances(const point& p, double& near_edge_distance,
    double& far_edge_distance) {
  // First, find the nearest vertex.
  S2Cell cell = get_cell();
  near_edge_distance = 100.0;
  far_edge_distance = -100.0;
  for (int k = 0; k < 4; k++) {
    double costheta = point::s2point_to_point(cell.GetVertexRaw(k)).dot(p);
    if (1.0 - costheta * costheta < near_edge_distance) {
      near_edge_distance = 1.0 - costheta * costheta;
    }
    if (1.0 - costheta * costheta > far_edge_distance) {
      far_edge_distance = 1.0 - costheta * costheta;
    }
  }

  // Edges are defined by great circles.  If a point is between the any opposed
  // pair of great circles, then there's a good chance than the nearest point
  // is on one of those edges (otherwise, the nearest point is a vertex).
  // Finding this point is more expensive than finding the nearest vertex, so
  // we check containment first.
  point_vector edges;
  circle_ptr_vector edge_caps;
  edges.reserve(4);
  edge_caps.reserve(4);
  for (int k = 0; k < 4; k++) {
    point edge = point::s2point_to_point(cell.GetEdgeRaw(k));
    edges.push_back(edge);
    edge_caps.push_back(circle_bound::from_height(edge, 0.0));
  }

  // Now iterate over the edge pairs, check containment and calculate the
  // distance if contained.
  bool nearest_point_on_edge = false;
  for (int k = 0; k < 2; k++) {
    if ((edge_caps[k])->contains(p) && (edge_caps[k + 2])->contains(p)) {
      // To find the nearest point on an edge (x), we need to find the normal
      // to the edge great circle that runs through p.  To find this, we solve
      // the following equations:
      //
      // edge(k).dot(p.cross(x)) == 0 && edge(k).dot(x) == 0.
      //
      // (i.e., the great circle that runs through p and x needs to be normal
      // to the edge great circle and x needs to be on the edge great circle).
      // The solution is x = edge(k).cross(p.cross(edge(k))).
      // point nearest_point = edges[k].cross(p.cross(edges[k])); // normalize?
      point nearest_point = edges[k].cross(point::cross(p, edges[k]));
      double costheta = nearest_point.dot(p);
      if (1.0 - costheta * costheta < near_edge_distance) {
        near_edge_distance = 1.0 - costheta * costheta;
        nearest_point_on_edge = true;
      }
      if (1.0 - costheta * costheta > far_edge_distance) {
        far_edge_distance = 1.0 - costheta * costheta;
      }

      nearest_point = edges[k + 2].cross(point::cross(p, edges[k + 2]));
      // nearest_point = edges[k + 2].cross(p.cross(edges[k + 2]));
      // normalize?
      costheta = nearest_point.dot(p);
      if (1.0 - costheta * costheta < near_edge_distance) {
        near_edge_distance = 1.0 - costheta * costheta;
        nearest_point_on_edge = true;
      }
      if (1.0 - costheta * costheta > far_edge_distance) {
        far_edge_distance = 1.0 - costheta * costheta;
      }
    }
  }

  return nearest_point_on_edge;
}

double pixel::nearest_edge_distance(const point& p) {
  double sin2theta_min, sin2theta_max;
  edge_distances(p, sin2theta_min, sin2theta_max);

  return sin2theta_min;
}

double pixel::farthest_edge_distance(const point& p) {
  double sin2theta_min, sin2theta_max;
  edge_distances(p, sin2theta_min, sin2theta_max);

  return sin2theta_max;
}

void pixel::neighbors(pixel_vector* pixels) const {
  neighbors(level(), pixels);
}

void pixel::neighbors(int level, pixel_vector* pixels) const {
  if (level > MAX_LEVEL || level < 0) {
    return;
  }

  vector<S2CellId> neighbor_ids;
  id_.AppendAllNeighbors(level, &neighbor_ids);
  pixels->clear();
  pixels->reserve(neighbor_ids.size());
  for (unsigned int k = 0; k < neighbor_ids.size(); k++) {
    pixels->push_back(pixel(neighbor_ids[k]));
  }
}

double pixel::contained_area(const pixel& pix) const {
  if (contains(pix)) {
    return pix.exact_area();
  } else if (pix.contains(*this)) {
    return exact_area();
  }
  return 0.0;
}

point pixel::get_center() const {
  return center_point();
}

circle_bound pixel::get_bound() const {
  circle_bound bound = circle_bound(center_point(), 0.0);
  for (int k = 0; k < 4; k++) {
    bound.add_point(vertex(k));
  }
  return bound;
}

void pixel::get_covering(pixel_vector* pixels) const {
  if (!pixels->empty())
    pixels->clear();
  pixels->push_back(*this);
}

void pixel::get_size_covering(const long max_pixels, pixel_vector* pixels) const {
  get_covering(pixels);
}

void pixel::get_area_covering(double fractional_area_tolerance,
    pixel_vector* pixels) const {
  get_covering(pixels);
}

void pixel::get_interior_covering(int max_level, pixel_vector* pixels) const {
  if (!pixels->empty())
    pixels->clear();
  if (level() <= max_level) {
    get_covering(pixels);
  }
}

void pixel::get_simple_covering(int level, pixel_vector* pixels) const {
  if (!pixels->empty())
    pixels->clear();

  if (level == this->level()) {
    pixels->push_back(*this);
  } else if (level > this->level()) {
    children(level, pixels);
  } else {
    pixels->push_back(parent(level));
  }
}

void pixel::get_center_covering(int level, pixel_vector* pixels) const {
  get_simple_covering(level, pixels);
}

} // end namespace s2omp
