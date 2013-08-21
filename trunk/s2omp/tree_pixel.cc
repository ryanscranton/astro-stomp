/*
 * tree_pixel.cc
 *
 *  Created on: Jul 8, 2012
 *      Author: cbmorrison
 */

// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)
// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains a variant on the Pixel class.  The goal here is to
// use the hierarchical nature of the Pixel class to form the basis for a
// spatial quad tree structure.  Hence, a given TreePixel object will have a
// number of points associated with it, where they have been stored in such a
// way that pair finding and K nearest neighbor searches will run in ln(N) time.
#include "core.h"
#include "tree_pixel.h"
#include "annulus_bound.h"

namespace s2omp {

tree_pixel::tree_pixel() {
  // The default constructor is an invalid pixel.
  set_id(0);
  initialize_node(DEFAULT_MAX_POINTS);
}

tree_pixel::tree_pixel(uint64 id) {
  set_id(id);
  initialize_node(DEFAULT_MAX_POINTS);
}

tree_pixel::tree_pixel(uint64 id, uint max_points) {
  set_id(id);
  initialize_node(max_points);
}

tree_pixel::~tree_pixel() {
  clear();
}

void tree_pixel::initialize_node(uint max_points) {
  weight_ = 0.0;
  maximum_points_ = 0;
  point_count_ = 0;
  subnodes_.clear();
  points_.clear();
}

tree_pixel* tree_pixel::from_point(const point& p, int level, uint max_points) {
  return new tree_pixel(point::point_to_id(p, level), max_points);
}

tree_pixel* tree_pixel::from_pixel(const pixel& pix, uint max_points) {
  return new tree_pixel(pix.id(), max_points);
}

// Moving private method here since it's referenced in add_point.
bool tree_pixel::initialize_subnodes() {
  if (is_leaf()) {
    return false;
  }

  subnodes_.reserve(4);
  for (pixel c = child_begin(); c != child_end(); c = c.next()) {
    subnodes_.push_back(from_pixel(c, maximum_points_));
  }

  for (point_ptr_iterator p_iter = points_.begin(); p_iter != points_.end(); ++p_iter) {
    bool transferred_point = false;
    for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); ++iter) {
      if ((*iter)->add_point(*p_iter)) {
        transferred_point = true;
        break;
      }
    }
    if (!transferred_point) {
      std::cout << "s2omp::tree_pixel::initialize_subnodes - "
          << "Failed to transfer point to any subnode.  Exiting.\n";
      exit(2);
    }
  }
  points_.clear();

  return !subnodes_.empty() && points_.empty();
}

bool tree_pixel::add_point(point* p) {
  // If the point is outside the pixel, ignore it.
  if (!contains(*p)) {
    return false;
  }

  // Add the point to the pixel if we are below capacity or if we're a leaf
  // node.
  if (point_count_ < maximum_points_ || is_leaf()) {
    if (point_count_ == 0)
      points_.reserve(maximum_points_);
    points_.push_back(p);
    weight_ += p->weight();
    point_count_++;

    return true;
  } else {
    // If we're at capacity, then this point will be added to a subnode.  Before
    // adding to any subnodes, we need to make sure that we've initialized them.
    if (!has_nodes()) {
      if (!initialize_subnodes()) {
        std::cout << "s2omp::tree_pixel::add_point - "
            << "Failed to initialize subnodes.  Exiting.\n";
        exit(2);
      }
    }

    // Iterate through the subnodes and add the point to the one that contains
    // it.
    for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); iter++) {
      if ((*iter)->contains(*p)) {
        // Provided that we successully
        if ((*iter)->add_point(p)) {
          weight_ += p->weight();
          point_count_++;

          return true;
        }
      }
    }
  }

  // If we've reached this point, then somehow we've failed to add the point
  // to either this node or any of the subnodes.
  return false;
}

bool tree_pixel::add_point(const point& p) {
  return add_point(point::copy_point(p));
}

// Moving these private methods here since the following methods use them.
void tree_pixel::direct_pair_count(const annulus_bound& bound,
    pair_weight* pairs) const {
  for (point_ptr_iterator iter = points_.begin(); iter != points_.end(); ++iter) {
    if (bound.contains(*(*iter))) {
      pairs->n_pairs++;
      pairs->total_weight += (*iter)->weight();
    }
  }
}

void tree_pixel::_find_pairs_recursion(const annulus_bound& bound,
    pair_weight* pairs) const {
  // Four possible cases:
  //   * If the bound doesn't intersect this node, we have no pairs.
  if (!bound.may_intersect(*this)) {
    return;
  }

  //   * If the bound contains the node, then every point in the node is a pair.
  if (bound.contains(*this)) {
    pairs->n_pairs += point_count_;
    pairs->total_weight += weight_;
    return;
  }

  //   * If we have no subnodes, then return a direct count of all the pairs.
  if (!points_.empty()) {
    direct_pair_count(bound, pairs);
    return;
  }

  //   * If we have subnodes, then recurse the bound through them and aggregate.
  for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); iter++) {
    (*iter)->_find_pairs_recursion(bound, pairs);
  }
}

long tree_pixel::find_pairs(const annulus_bound& bound) const {
  pair_weight pairs;
  _find_pairs_recursion(bound, &pairs);

  return pairs.n_pairs;
}

double tree_pixel::find_weighted_pairs(const annulus_bound& bound) const {
  pair_weight pairs;
  _find_pairs_recursion(bound, &pairs);

  return pairs.total_weight;
}

void tree_pixel::_neighbor_recursion(const point& p, tree_neighbor* neighbors) const {
  neighbors->add_node();

  if (!points_.empty()) {
    // We have no subnodes in this tree, so we'll just iterate over the
    // points here and take the nearest N neighbors.
    for (point_ptr_iterator iter = points_.begin(); iter != points_.end(); ++iter) {
      neighbors->test_point(*iter);
    }
    return;
  }

  // This node is the root node for our tree, so we first find the sub-node
  // that contains the point and start recursing there.
  //
  // While we iterate through the nodes, we'll also calculate the edge
  // distances for those nodes that don't contain the point and store them
  // in a priority queue.  This will let us do a follow-up check on nodes in
  // the most productive order.
  pixel_queue pix_queue;
  for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); ++iter) {
    if ((*iter)->contains(p)) {
      (*iter)->_neighbor_recursion(p, neighbors);
    } else {
      distance_pixel_pair dist_pair((*iter)->nearest_edge_distance(p), (*iter));
      pix_queue.push(dist_pair);
    }
  }

  // That should give us back a tree_neighbor object that contains a workable
  // set of neighbors and a search radius for possible matches.  Now we just
  // need to iterate over those sub-nodes that didn't contain the input point
  // to verify that there can't be any points in their sub-nodes which might
  // be closer to the input point.
  //
  // There's also the possibility that the input point is completely outside
  // our tree.  In that case (where the number of neighbors in the
  // TreeNeighbor object is less than the maximum), we want to check
  // all nodes.
  while (!pix_queue.empty()) {
    double pix_distance = pix_queue.top().first;
    tree_pixel* pix_iter = pix_queue.top().second;
    if (pix_distance < neighbors->max_distance()) {
      pix_iter->_neighbor_recursion(p, neighbors);
    }
    pix_queue.pop();
  }
}

long tree_pixel::find_k_nearest_neighbors(const point& p, uint n_neighbors,
    point_vector* neighbor_vector) const {
  tree_neighbor neighbors(p, n_neighbors);

  _neighbor_recursion(p, &neighbors);

  neighbors.nearest_neighbors(neighbor_vector, false);

  return neighbors.nodes_visited();
}

long tree_pixel::find_nearest_neighbor(const point& p, point* neighbor) const {
  point_vector p_vector;

  long nodes_visited = find_k_nearest_neighbors(p, 1, &p_vector);

  neighbor = point::copy_point(p_vector[0]);

  return nodes_visited;
}

double tree_pixel::k_nearest_neighbor_distance(const point& p,
    uint n_neighbors, long& nodes_visited) const {

  tree_neighbor neighbors(p, n_neighbors);

  _neighbor_recursion(p, &neighbors);

  nodes_visited = neighbors.nodes_visited();

  return neighbors.max_angular_distance();
}

double tree_pixel::nearest_neighbor_distance(const point& p,
    long& nodes_visited) const {
  return k_nearest_neighbor_distance(p, 1, nodes_visited);
}

bool tree_pixel::closest_match(const point& p, double max_angular_distance,
    point* match) const {

  tree_neighbor neighbors(p, 1, max_angular_distance);

  _neighbor_recursion(p, &neighbors);

  bool found_match = false;
  if (neighbors.n_neighbors() == neighbors.max_neighbors()
      && neighbors.max_angular_distance() < max_angular_distance) {
    found_match = true;

    point_vector p_vector;
    neighbors.nearest_neighbors(&p_vector, false);
    match = point::copy_point(p_vector[0]);
  }

  return found_match;
}

long tree_pixel::n_points(const pixel& pix) const {
  // If the input pixel contains this node, then all points are contained.
  if (pix.contains(*this)) {
    return point_count_;
  }

  // If that's not the case and this node doesn't contain the input pixel,
  // then no points are contained.
  if (!contains(pix)) {
    return 0;
  }

  // Finally, if the input pixel is a child of this node, then we need to
  // determine how many points in this node are contained.  Start with the case
  // where we have points in this node.
  long contained_points = 0;
  if (!points_.empty()) {
    for (point_ptr_iterator iter = points_.begin(); iter != points_.end(); iter++) {
      if (pix.contains(*(*iter)))
        contained_points++;
    }

    return contained_points;
  }

  // If we have subnodes, iterate over them and return the aggregate.
  for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); iter++) {
    contained_points += (*iter)->n_points(pix);
  }

  return contained_points;
}

double tree_pixel::weight(const pixel& pix) const {
  // If the input pixel contains this node, then all points are contained.
  if (pix.contains(*this)) {
    return weight_;
  }

  // If that's not the case and this node doesn't contain the input pixel, then
  // no points are contained.
  if (!contains(pix)) {
    return 0.0;
  }

  // Finally, if the input pixel is a child of this node, then we need to
  // determine how many points in this node are contained.  Start with the case
  // where we have points in this node.
  double total_weight = 0.0;
  if (!points_.empty()) {
    for (point_ptr_iterator iter = points_.begin(); iter != points_.end(); iter++) {
      if (pix.contains(*(*iter)))
        total_weight = (*iter)->weight();
    }

    return total_weight;
  }

  // If we have subnodes, iterate over them and return the aggregate.
  for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); iter++) {
    total_weight += (*iter)->weight(pix);
  }

  return total_weight;
}

double tree_pixel::covering_fraction() const {
  // If the node is empty, then the covering fraction is 0.
  if (point_count_ == 0) {
    return 0.0;
  }

  // If the node contains only points and no subnodes, then the covering
  // fraction is unity.
  if (!points_.empty()) {
    return 1.0;
  }

  // Otherwise, we recursively count the number of subnodes in the tree that
  // contain points.
  double total_area = exact_area();
  double covered_area = 0.0;

  for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); iter++) {
    covered_area += (*iter)->covering_fraction() * (*iter)->exact_area();
  }

  return covered_area / total_area;
}

double tree_pixel::covering_fraction(const pixel& pix) const {
  // If the input pixel contains this node, then the covering fraction is the
  // covering_fraction of this node scaled by the relative areas.
  if (pix.contains(*this)) {
    return exact_area() * covering_fraction() / pix.exact_area();
  }

  // If that's not true and this node does not contain the input pixel, then the
  // covering fraction is 0.  This also holds if this node is empty.
  if (!contains(pix) || point_count_ == 0) {
    return 0.0;
  }

  // Finally, we deal with the case that our node contains either points or
  // subnodes and the input pixel.  First, deal with the former case.  If any
  // of our contained points is also contained by the input pixel, then the
  // contained fraction is unity.  Otherwise, the fraction is zero.
  if (!points_.empty()) {
    for (point_ptr_iterator iter = points_.begin(); iter != points_.end(); iter++) {
      if (pix.contains(*(*iter))) {
        return 1.0;
      }
    }

    return 0.0;
  }

  // If we have subnodes, then iterate over them and return the aggregate.
  double covered_fraction = 0.0;
  for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); iter++) {
    covered_fraction += (*iter)->covering_fraction(pix);
  }

  return covered_fraction;
}

void tree_pixel::points(point_vector* p) const {
  if (!p->empty())
    p->clear();

  if (point_count_ == 0) {
    return;
  }

  // The output vector should be the same size as our current point count.
  p->reserve(point_count_);

  // If we have points in this node, then copy them to the output vector.
  if (!points_.empty()) {
    for (point_ptr_iterator iter = points_.begin(); iter != points_.end(); iter++) {
      p->push_back(point::copy_point(*iter));
    }
    return;
  }

  // If we contain subnodes, then we need to copy them from the subnodes
  // to the output vector.  This should probably be avoid if at all possible.
  for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); iter++) {
    point_vector tmp_points;
    (*iter)->points(&tmp_points);
    for (point_iterator p_iter = tmp_points.begin(); p_iter != tmp_points.end(); p_iter++) {
      p->push_back(*p_iter);
    }
  }
}

void tree_pixel::points(const pixel& pix, point_vector* p) const {
  if (!p->empty())
    p->clear();

  // If the input pixel contains this node, then copy all of this node's points
  // into the output array.
  if (pix.contains(*this)) {
    points(p);
    return;
  }

  // If that's not true and this node does not contain the input pixel, then
  // there are no points to copy.  This also holds if this node is empty.
  if (!contains(pix) || point_count_ == 0) {
    return;
  }

  // If the input pixel is contained in the node, then we need to copy the
  // contained points to the output array.  Start with the case where our
  // node contains points and copy over the ones contained in the input pixel.
  p->reserve(point_count_);

  if (!points_.empty()) {
    for (point_ptr_iterator iter = points_.begin(); iter != points_.end(); iter++) {
      if (pix.contains(*(*iter))) {
        p->push_back(point::copy_point(*iter));
      }
    }
    return;
  }

  // If we contain subnodes, then we need to copy them from the subnodes
  // to the output vector.  This should probably be avoid if at all possible.
  for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); iter++) {
    point_vector tmp_points;
    (*iter)->points(pix, &tmp_points);
    for (point_iterator p_iter = tmp_points.begin(); p_iter != tmp_points.end(); p_iter++) {
      p->push_back(*p_iter);
    }
  }
}

long tree_pixel::n_nodes() const {
  if (!points_.empty()) {
    return 1;
  }

  long total_nodes = 1;
  for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); iter++) {
    total_nodes += (*iter)->n_nodes();
  }

  return total_nodes;
}

void tree_pixel::clear() {
  for (point_ptr_iterator iter = points_.begin(); iter != points_.end(); iter++) {
    delete *iter;
  }
  points_.clear();

  for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); ++iter) {
    (*iter)->clear();
    delete *iter;
  }
  subnodes_.clear();
}

tree_neighbor::tree_neighbor(const point& reference_point) {
  reference_point_ = reference_point;
  n_neighbors_ = DEFAULT_N_NEIGHBORS;
  max_distance_ = sin(DEG_TO_RAD * DEFAULT_MAX_NEIGHBOR_DISTANCE);
  max_distance_ *= max_distance_;
  n_nodes_visited_ = 0;
}

tree_neighbor::tree_neighbor(const point& reference_point, uint n_neighbors) {
  reference_point_ = reference_point;
  n_neighbors_ = n_neighbors;
  max_distance_ = sin(DEG_TO_RAD * DEFAULT_MAX_NEIGHBOR_DISTANCE);
  max_distance_ *= max_distance_;
  n_nodes_visited_ = 0;
}

tree_neighbor::tree_neighbor(const point& reference_point, uint n_neighbors,
    double max_angular_distance) {
  reference_point_ = reference_point;
  n_neighbors_ = n_neighbors;
  max_distance_ = sin(DEG_TO_RAD * max_angular_distance);
  max_distance_ *= max_distance_;
  n_nodes_visited_ = 0;
}

void tree_neighbor::nearest_neighbors(point_vector* p, bool save_neighbors) {
  if (!p->empty())
    p->clear();
  std::vector<distance_point_pair> backup_copy;

  while (!point_queue_.empty()) {
    distance_point_pair dist_pair = point_queue_.top();
    point_queue_.pop();

    p->push_back(point::copy_point(dist_pair.second));
    backup_copy.push_back(dist_pair);
  }

  if (save_neighbors) {
    for (uint i = 0; i < backup_copy.size(); i++) {
      point_queue_.push(backup_copy[i]);
    }
  }
}

point tree_neighbor::nearest_neighbor() const {
  return point::copy_point(point_queue_.top().second);
}

bool tree_neighbor::test_point(point* test_point) {
  double costheta = reference_point_.dot(*test_point);
  double sin2theta = 1.0 - costheta * costheta;

  if (sin2theta < max_distance_ || n_neighbors() < max_neighbors()) {
    if (n_neighbors() == max_neighbors())
      point_queue_.pop();
    point_queue_.push(distance_point_pair(sin2theta, test_point));
    max_distance_ = point_queue_.top().first;

    return true;
  }

  return false;
}

double tree_neighbor::max_angular_distance() {
  return RAD_TO_DEG * asin(sqrt(fabs(max_distance_)));
}

} //end namespace s2omp

