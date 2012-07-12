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

#include "core.h"
#include "pixel.h"
#include "point.h"
#include "circle_bound.h"

pixel* pixel::from_point(point& p) {
  return new pixel(p.id());
}

static pixel* from_point(const point& p, int level) {
  return new pixel(p.id(level));
}

pixel pixel::parent() const {
  return pixel(id_.parent().id());
}

pixel pixel::parent(int level) const {
  return pixel(id_.parent(level).id());
}

void pixel::children(pixel_vector* child_pixels) const {
  child_pixels->clear();
  if (!is_leaf()) {
    child_pixels->reserve(4);
    for (int k = 0; k < 4; k++) {
      child_pixels->push_back(pixel(id_.child(k).id()));
    }
  }
}

void pixel::children(int level, pixel_vector* child_pixels) const {
  child_pixels->clear();
  if (!is_leaf() && level <= MAX_LEVEL) {
    for (pixel c = pix.child_begin(level); c != pix.child_end(level); c
        = c.next()) {
      child_pixels->push_back(c);
    }
  }
}

pixel pixel::child_begin() const {
  return pixel(id_.child_begin().id());
}

pixel pixel::child_begin(int level) const {
  return pixel(id_.child_begin(level).id());
}

pixel pixel::child_end() const {
  return pixel(id_.child_end().id());
}

pixel pixel::child_end(int level) const {
  return pixel(id_.child_end(level).id());
}

pixel pixel::next() const {
  return pixel(id_.next().id());
}

pixel pixel::prev() const {
  return pixel(id_.prev().id());
}

pixel pixel::next_wrap() const {
  return pixel(id_.next_wrap().id());
}

pixel pixel::prev_wrap() const {
  return pixel(id_.prev_wrap().id());
}

S2::S2Cell pixel::get_cell() const {
  return S2::S2Cell(id_);
}

point pixel::s2point_to_point(const S2::S2Point& p) const {
  return point(p.x(), p.y(), p.z(), 1.0);
}

S2::S2Point pixel::point_to_s2point(const point& p) const {
  return S2::S2Point(p.unit_sphere_x(), p.unit_sphere_y(), p.unit_sphere_z());
}

double pixel::average_area(int level) {
  return S2::S2Cell.AverageArea(level) * STRAD_TO_DEG2;
}

double pixel::average_area() const {
  return get_cell().AverageArea() * STRAD_TO_DEG2;
}

double pixel::exact_area() const {
  return get_cell().ApproxArea() * STRAD_TO_DEG2;
}

bool pixel::contains(const point& p) const {
  return get_cell().Contains(point_to_s2point(p));
}

bool pixel::contains(const pixel& pix) const {
  return get_cell().Contains(S2CellId(pix.id()));
}

point pixel::center_point() const {
  return s2point_to_point(id_.ToPointRaw()); // point will normalize itself.
}

point pixel::vertex(int k) const {
  return s2point_to_point(get_cell().GetVertexRaw(k));
}

point pixel::edge(int k) const {
  return s2point_to_point(get_cell().GetEdgeRaw(k));
}

bool pixel::edge_distances(point& p, double& near_edge_distance,
    double& far_edge_distance) {
  // First, find the nearest vertex.
  S2::S2Cell cell = get_cell();
  double near_edge_distance = 100.0;
  double far_edge_distance = -100.0;
  for (int k = 0; k < 4; k++) {
    double costheta = s2point_to_point(cell.GetVertexRaw(k)).dot(p);
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
  circle_vector edge_caps;
  edges.reserve(4);
  edge_caps.reserve(4);
  for (int k = 0; k < 4; k++) {
    point edge = s2point_to_point(cell.getRawEdge(k));
    edges.push_back(edge);
    edge_caps.push_back(circle_bound.from_height(edge, 0.0));
  }

  // Now iterate over the edge pairs, check containment and calculate the
  // distance if contained.
  bool nearest_point_on_edge = false;
  for (int k = 0; k < 2; k++) {
    if (edge_caps[k].contains(p) && edge_caps[k + 2].contains(p)) {
      // To find the nearest point on an edge, we need to find the normal to
      // the edge great circle that runs through p.  To find this, we solve
      // the following equations:
      //
      // edge(k).dot(p.cross(x)) == 0 && edge(k).dot(x) == 0.
      //
      // The solution is x = edge(k).cross(p.cross(edge(k))).
      point nearest_point = edge[k].cross(p.cross(edge[k])); // normalize?
      double costheta = nearest_point.dot(p);
      if (1.0 - costheta * costheta < near_edge_distance) {
        near_edge_distance = 1.0 - costheta * costheta;
        nearest_point_on_edge = true;
      }
      if (1.0 - costheta * costheta > far_edge_distance) {
        far_edge_distance = 1.0 - costheta * costheta;
      }

      nearest_point = edge[k + 2].cross(p.cross(edge[k + 2])); // normalize?
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

double pixel::nearest_edge_distance(point& p) {
  double sin2theta_min, sin2theta_max;
  edge_distances(p, sin2theta_min, sin2theta_max);

  return sin2theta_min;
}

double pixel::farthest_edge_distance(point& p) {
  double sin2theta_min, sin2theta_max;
  edge_distances(p, sin2theta_min, sin2theta_max);

  return sin2theta_max;
}

void pixel::neighbors(pixel_vector* pixels) const {
  neighbors(level(), pixels);
}

void pixel::neighbors(int level, pixel_vector* pixels) const {
  if (level <= MAX_LEVEL) {
    vector<S2::S2CellId> neighbor_ids;
    id_.AppendAllNeighbors(level, &neighbor_ids);
    pixels->clear();
    for (int k = 0; k < neighbor_ids.size(); k++) {
      pixels->push_back(pixel(neighbor_ids[k].id()));
    }
  }
}

point pixel::get_random_point() const {
  return point(0.0, 0.0, 1.0, 1.0);
}

void pixel::get_random_points(long n_points, pixel_vector* points) const {
  return;
}
