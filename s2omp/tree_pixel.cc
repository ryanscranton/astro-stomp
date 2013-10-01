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
  initialize_node(DEFAULT_NODE_CAPACITY);
  cell_ = S2Cell(S2CellId(id()));
}

tree_pixel::tree_pixel(uint64 id) {
  set_id(id);
  initialize_node(DEFAULT_NODE_CAPACITY);
  cell_ = S2Cell(S2CellId(id));
}

tree_pixel::tree_pixel(uint64 id, uint node_capacity) {
  set_id(id);
  initialize_node(node_capacity);
  cell_ = S2Cell(S2CellId(id));
}

tree_pixel::~tree_pixel() {
  clear();
}

void tree_pixel::initialize_node(uint node_capacity) {
  node_capacity_ = node_capacity;
  weight_ = 0.0;
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
    subnodes_.push_back(from_pixel(c, node_capacity_));
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
  if (point_count_ < node_capacity_ || is_leaf()) {
    if (point_count_ == 0)
      points_.reserve(node_capacity_);
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
  for (int k = 0; k < subnodes_.size(); k++) {
    tree_pixel* node = subnodes_[k];
    if (node->contains(p)) {
      node->_neighbor_recursion(p, neighbors);
    } else {
      double dist = node->nearest_edge_distance(p);
      if (dist < neighbors->max_distance() || neighbors->n_neighbors() == 0) {
        pix_queue.push(distance_pixel_pair(dist, node));
      }
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
  // tree_neighbor object is less than the maximum), we want to check
  // all nodes, starting with the node nearest the input point.
  while (!pix_queue.empty()) {
    double pix_distance = pix_queue.top().first;
    tree_pixel* node = pix_queue.top().second;
    if (pix_distance < neighbors->max_distance()) {
      node->_neighbor_recursion(p, neighbors);
    }
    pix_queue.pop();
  }
}

void tree_pixel::find_k_nearest_neighbors(const point& p, uint n_neighbors,
    point_vector* neighbor_vector) const {
  tree_neighbor neighbors(p, n_neighbors);

  _neighbor_recursion(p, &neighbors);

  neighbors.nearest_neighbors(neighbor_vector, false);
}

point tree_pixel::find_nearest_neighbor(const point& p) const {
  point_vector points;

  find_k_nearest_neighbors(p, 1, &points);

  return points[0];
}

double tree_pixel::k_nearest_neighbor_distance(const point& p,
    uint n_neighbors) const {

  tree_neighbor neighbors(p, n_neighbors);

  _neighbor_recursion(p, &neighbors);

  return neighbors.max_distance_deg();
}

double tree_pixel::nearest_neighbor_distance(const point& p) const {
  return k_nearest_neighbor_distance(p, 1);
}

bool tree_pixel::closest_match(const point& p, double max_distance_deg,
    point& match) const {

  tree_neighbor neighbors(p, 1, max_distance_deg);

  _neighbor_recursion(p, &neighbors);

  if (neighbors.n_neighbors() != neighbors.max_neighbors()
      || neighbors.max_distance_deg() > max_distance_deg) {
    return false;
  }

  match = neighbors.nearest_neighbor();

  return true;
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
  double covered_area = 0.0;

  for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); iter++) {
    covered_area += (*iter)->covering_fraction() * (*iter)->exact_area();
  }

  return covered_area / exact_area();
}

double tree_pixel::covering_fraction(const pixel& pix) const {
  // If the input pixel contains this node, then the covering fraction is the
  // covering_fraction of this node scaled by the relative areas.
  if (pix.contains(*this)) {
    return exact_area() * covering_fraction() / pix.exact_area();
  }

  // If that's not true and this node does not contain the input pixel, then the
  // covering fraction is 0.  This also holds if this node is empty.
  if (!contains(pix) || is_empty()) {
    return 0.0;
  }

  // Finally, we deal with the case that our node contains either points or
  // subnodes and the input pixel.  First, deal with the latter case.  If we
  // have subnodes, then iterate over them and return the value for the
  // containing subnode.
  if (!subnodes_.empty()) {
    for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); iter++) {
      if ((*iter)->contains(pix)) {
        return (*iter)->covering_fraction(pix);
      }
    }
  }

  // If any of our contained points is also contained by the input pixel, then
  // the contained fraction is unity.  Otherwise, the fraction is zero.
  for (point_ptr_iterator iter = points_.begin(); iter != points_.end(); iter++) {
    if (pix.contains(*(*iter))) {
      return 1.0;
    }
  }

  return 0.0;
}

void tree_pixel::copy_points(point_vector* points) const {
  if (!points->empty())
    points->clear();

  if (point_count_ == 0) {
    return;
  }

  // The output vector should be the same size as our current point count.
  points->reserve(point_count_);

  // If we have points in this node, then copy them to the output vector.
  if (!points_.empty()) {
    for (point_ptr_iterator iter = points_.begin(); iter != points_.end(); iter++) {
      points->push_back(point::copy_point(*iter));
    }
    return;
  }

  // If we contain subnodes, then we need to copy them from the subnodes
  // to the output vector.  This should probably be avoid if at all possible.
  for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); iter++) {
    point_vector tmp_points;
    (*iter)->copy_points(&tmp_points);
    points->insert(points->end(), tmp_points.begin(), tmp_points.end());
  }
}

void tree_pixel::copy_points(const pixel& pix, point_vector* points) const {
  if (!points->empty())
    points->clear();

  // If the input pixel contains this node, then copy all of this node's points
  // into the output array.
  if (pix.contains(*this)) {
    copy_points(points);
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
  points->reserve(point_count_);

  if (!points_.empty()) {
    for (point_ptr_iterator iter = points_.begin(); iter != points_.end(); iter++) {
      if (pix.contains(*(*iter))) {
        points->push_back(point::copy_point(*iter));
      }
    }
    return;
  }

  // If we contain subnodes, then we need to copy them from the subnodes
  // to the output vector.  This should probably be avoid if at all possible.
  for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); iter++) {
    point_vector tmp_points;
    (*iter)->copy_points(pix, &tmp_points);
    points->insert(points->end(), tmp_points.begin(), tmp_points.end());
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

  for (tree_ptr_iterator iter = subnodes_.begin(); iter != subnodes_.end(); ++iter) {
    (*iter)->clear();
    delete *iter;
  }

  initialize_node(node_capacity_);
}

S2Cell tree_pixel::get_cell() const {
  return cell_;
}

tree_neighbor::tree_neighbor(const point& reference_point) {
  reference_point_ = reference_point;
  max_neighbors_ = DEFAULT_MAX_NEIGHBORS;
  max_distance_ = sin(DEG_TO_RAD * DEFAULT_MAX_NEIGHBOR_DISTANCE);
  max_distance_ *= max_distance_;
  n_nodes_visited_ = 0;
}

tree_neighbor::tree_neighbor(const point& reference_point, uint max_neighbors) {
  reference_point_ = reference_point;
  max_neighbors_ = max_neighbors;
  max_distance_ = sin(DEG_TO_RAD * DEFAULT_MAX_NEIGHBOR_DISTANCE);
  max_distance_ *= max_distance_;
  n_nodes_visited_ = 0;
}

tree_neighbor::tree_neighbor(const point& reference_point, uint max_neighbors,
    double max_distance_deg) {
  reference_point_ = reference_point;
  max_neighbors_ = max_neighbors;
  max_distance_ = sin(DEG_TO_RAD * max_distance_deg);
  max_distance_ *= max_distance_;
  n_nodes_visited_ = 0;
}

tree_neighbor::~tree_neighbor() {
  while (!point_queue_.empty()) point_queue_.pop();
}

void tree_neighbor::nearest_neighbors(point_vector* points,
    bool save_neighbors) {
  if (!points->empty()) points->clear();
  std::vector<distance_point_pair> backup_copy;

  points->reserve(point_queue_.size());
  while (!point_queue_.empty()) {
    distance_point_pair dist_pair = point_queue_.top();
    point_queue_.pop();

    points->push_back(point::copy_point(dist_pair.second));
    backup_copy.push_back(dist_pair);
  }

  if (save_neighbors) {
    for (uint i = 0; i < backup_copy.size(); i++) {
      point_queue_.push(backup_copy[i]);
    }
  }
}

point tree_neighbor::nearest_neighbor() {
  point_vector points;
  nearest_neighbors(&points, true);
  return points.back();
}

bool tree_neighbor::test_point(point* test_point) {
  double sin2theta = reference_point_.cross_norm2(*test_point);

  if (sin2theta < max_distance_ || n_neighbors() < max_neighbors()) {
    if (n_neighbors() == max_neighbors()) {
      point_queue_.pop();
    }
    point_queue_.push(distance_point_pair(sin2theta, test_point));
    max_distance_ = point_queue_.top().first;

    return true;
  }

  return false;
}

double tree_neighbor::max_distance_deg() {
  return RAD_TO_DEG * sqrt(asin(fabs(max_distance_)));
}

} //end namespace s2omp

