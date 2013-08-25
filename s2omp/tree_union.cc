/*
 * tree_union.cc
 *
 *  Created on: Jul 16, 2012
 *      Author: cbmorrison
 */

#include "tree_union.h"
#include "angular_bin-inl.h"
#include "annulus_bound.h"
#include "point.h"
#include "region_map.h"

namespace s2omp {

tree_union::tree_union() {
  maximum_points_ = DEFAULT_MAX_POINTS;
  level_ = DEFAULT_LEVEL;
  point_count_ = 0;
  weight_ = 0.0;
  area_ = 0.0;
  modified_ = false;
  initialized_bound_ = false;
}

tree_union::tree_union(int level) {
  maximum_points_ = DEFAULT_MAX_POINTS;
  level_ = level;
  point_count_ = 0;
  weight_ = 0.0;
  area_ = 0.0;
  modified_ = false;
  initialized_bound_ = false;
}

tree_union::tree_union(int level, int max_points) {
  maximum_points_ = max_points;
  point_count_ = 0;
  level_ = level;
  weight_ = 0.0;
  area_ = 0.0;
  modified_ = false;
  initialized_bound_ = false;
}

// Moving private method here since add_point uses it.
tree_map_iterator tree_union::add_node(uint64 id) {
  tree_map_insert_iterator iter = tree_map_.insert(std::pair<uint64,
      tree_pixel*>(id, new tree_pixel(id, maximum_points_)));
  if (!iter.second) {
    std::cout << "s2omp::tree_union::add_node - "
        << "Creating new tree_union node failed. Exiting.\n";
    exit(2);
  }

  nodes_.insert(pixel(id));
  modified_ = true;

  return iter.first;
}

bool tree_union::add_point(const point& p) {
  tree_map_iterator iter = tree_map_.find(p.id(level_));

  if (iter == tree_map_.end()) {
    // Add a new node to contain the point.
    iter = add_node(p.id(level_));
  }

  if (!(*iter).second->add_point(p)) {
    std::cout << "s2omp::tree_union::add_point - "
        << "Adding point to tree_union failed. Exiting.\n";
    exit(2);
  }

  point_count_++;
  weight_ += p.weight();
  modified_ = true;

  return modified_;
}

long tree_union::find_pairs(const annulus_bound& bound) const {
  pixel_vector pixels;
  bound.get_simple_covering(level_, &pixels);

  long n_pairs = 0;
  for (pixel_iterator iter = pixels.begin(); iter != pixels.end(); ++iter) {
    tree_map_iterator tree_iter = tree_map_.find(iter->id());
    if (tree_iter == tree_map_.end())
      continue;

    n_pairs += (*tree_iter).second->find_pairs(bound);
  }

  return n_pairs;
}

double tree_union::find_weighted_pairs(const annulus_bound& bound) const {
  pixel_vector pixels;
  bound.get_simple_covering(level_, &pixels);

  double total_weight = 0.0;
  for (pixel_iterator iter = pixels.begin(); iter != pixels.end(); ++iter) {
    tree_map_iterator tree_iter = tree_map_.find(iter->id());
    if (tree_iter == tree_map_.end())
      continue;

    total_weight += (*tree_iter).second->find_weighted_pairs(bound);
  }

  return total_weight;
}

bool tree_union::find_pairs(const point_vector& p, angular_bin* bin) const {
  return find_pairs_with_regions(p, region_map(), bin);
}

bool tree_union::find_pairs_with_regions(const point_vector& points,
    const region_map& regions, angular_bin* bin) const {
  // If we're regionating the measurement, we need to make sure that the
  // regions are defined at a coarser resolution than our tree_union.
  if (!regions.is_empty() && regions.level() > level_) {
    std::cout << "s2omp::tree_union::find_pairs_with_regions - "
        << "Cannot use region map at " << regions.level() << " with level "
        << level_ << " tree_union\n";
    return false;;
  }

  pixel_vector pixels;
  int point_region = INVALID_REGION_VALUE;

  for (point_iterator point_iter = points.begin(); point_iter != points.end(); ++point_iter) {
    // Find the covering pixels for an annulus bound defined by this point
    // and the input angular_bin.
    annulus_bound* bound = annulus_bound::from_angular_bin(*point_iter, *bin);
    bound->get_simple_covering(level_, &pixels);

    // Set our point_region value based on the input region_map if we're using
    // it.  If we have a non-empty region_map and our region value is -1, then
    // we fail out of the iteration.
    point_region = region_map::find_region(regions, *point_iter);

    // Iterate over the covering pixels and find pairs for all of the
    for (pixel_iterator cover_iter = pixels.begin(); cover_iter != pixels.end(); ++cover_iter) {
      tree_map_iterator tree_iter = tree_map_.find(cover_iter->id());
      if (tree_iter == tree_map_.end())
        continue;

      pair_weight pairs;
      (*tree_iter).second->_find_pairs_recursion(*bound, &pairs);
      bin->add_to_pair_wtheta(pairs.total_weight, pairs.n_pairs, point_region,
          region_map::find_region(regions, *cover_iter));
    }

    delete bound;
  }

  return true;
}

void tree_union::_neighbor_recursion(const point& p, tree_neighbor* neighbors) const {
  // First we need to find out if the input point is within our map area.
  tree_map_iterator iter = tree_map_.find(p.id(level_));

  // If a node containing this point exists, then start finding neighbors there.
  if (iter != tree_map_.end())
    iter->second->_neighbor_recursion(p, neighbors);

  // That should give us back a tree_neighbor object that contains a workable
  // set of neighbors and a search radius for possible matches.  Now we just
  // need to iterate over those nodes that didn't contain the input point
  // to verify that there can't be any points in their sub-nodes which might
  // be closer to the input point.
  //
  // There's also the possibility that the input point is completely outside
  // our tree.  In that case (where the number of neighbors in the
  // TreeNeighbor object is less than the maximum), we want to check
  // all nodes.
  pixel_vector pixels;
  if (neighbors->n_neighbors() == neighbors->max_neighbors()) {
    // We've got a starting list of neighbors, so we only have to look at
    // nodes within our current bounding radius.
    circle_bound* bound = circle_bound::from_radius(p,
        neighbors->max_angular_distance());
    bound->get_simple_covering(level_, &pixels);
    delete bound;
  } else {
    for (pixel_set_iterator pix_iter = nodes_.begin(); pix_iter != nodes_.end(); ++pix_iter) {
      pixels.push_back(*pix_iter);
    }
  }

  // Now we construct a priority queue so that we're search the nodes closest
  // to the input point first.
  pixel_queue pix_queue;
  for (pixel_iterator pix_iter = pixels.begin(); pix_iter != pixels.end(); ++pix_iter) {
    iter = tree_map_.find(pix_iter->id());
    if (iter != tree_map_.end() && !pix_iter->contains(p)) {
      double min_edge_distance, max_edge_distance;
      iter->second->edge_distances(p, min_edge_distance, max_edge_distance);
      distance_pixel_pair dist_pair(min_edge_distance, iter->second);
      pix_queue.push(dist_pair);
    }
  }

  // And iterate over that queue to check for neighbors.
  while (!pix_queue.empty()) {
    double pix_distance = pix_queue.top().first;
    tree_pixel* pix_iter = pix_queue.top().second;
    if (pix_distance < neighbors->max_distance()) {
      pix_iter->_neighbor_recursion(p, neighbors);
    }
    pix_queue.pop();
  }
}

long tree_union::find_k_nearest_neighbors(const point& p, int n_neighbors,
    point_vector* points) const {
  tree_neighbor neighbors(p, n_neighbors);

  _neighbor_recursion(p, &neighbors);

  neighbors.nearest_neighbors(points, false);

  return neighbors.nodes_visited();
}

long tree_union::find_nearest_neighbor(const point& p, point* neighbor) const {
  point_vector p_vector;

  long nodes_visited = find_k_nearest_neighbors(p, 1, &p_vector);

  neighbor = point::copy_point(p_vector[0]);

  return nodes_visited;
}

double tree_union::k_nearest_neighbor_distance(const point& p, int n_neighbors,
    long& nodes_visited) const {
  tree_neighbor neighbors(p, n_neighbors);

  _neighbor_recursion(p, &neighbors);

  nodes_visited = neighbors.nodes_visited();

  return neighbors.max_angular_distance();
}

double tree_union::nearest_neighbor_distance(const point& p,
    long& nodes_visited) const {
  return k_nearest_neighbor_distance(p, 1, nodes_visited);
}

void tree_union::_match_recursion(const point& p, tree_neighbor* neighbors) const {
  // First we need to find out if the input point is within our map area.
  pixel center_pix = p.to_pixel(level_);
  tree_map_iterator iter = tree_map_.find(center_pix.id());

  // If a node containing this point exists, then start finding neighbors there.
  if (iter != tree_map_.end())
    iter->second->_neighbor_recursion(p, neighbors);

  // There's also a possibility that the matching point is just on the other
  // side of a pixel boundary.  To see if that's possible, check the edge
  // distance to our reference pixel.
  if (center_pix.nearest_edge_distance(p) > neighbors->max_distance()) {
    // We're far enough from any boundary that we don't need to look at any
    // other pixels.
    return;
  }

  pixel_vector pixels;
  circle_bound* bound = circle_bound::from_radius(p,
      neighbors->max_angular_distance());
  bound->get_simple_covering(level_, &pixels);
  delete bound;

  // Now we construct a priority queue so that we're search the nodes closest
  // to the input point first.
  pixel_queue pix_queue;
  for (pixel_iterator pix_iter = pixels.begin(); pix_iter != pixels.end(); ++pix_iter) {
    iter = tree_map_.find(pix_iter->id());
    if (iter != tree_map_.end() && !pix_iter->contains(p)) {
      double min_edge_distance, max_edge_distance;
      iter->second->edge_distances(p, min_edge_distance, max_edge_distance);
      distance_pixel_pair dist_pair(min_edge_distance, iter->second);
      pix_queue.push(dist_pair);
    }
  }

  // And iterate over that queue to check for neighbors.
  while (!pix_queue.empty()) {
    double pix_distance = pix_queue.top().first;
    tree_pixel* pix_iter = pix_queue.top().second;
    if (pix_distance < neighbors->max_distance()) {
      pix_iter->_neighbor_recursion(p, neighbors);
    }
    pix_queue.pop();
  }
}

bool tree_union::closest_match(const point& p, double max_angular_distance,
    point& match) const {
  tree_neighbor neighbors(p, 1, max_angular_distance);

  _match_recursion(p, &neighbors);

  // If we weren't able to find a nearest point or if the distance to the
  // neighbor point is greater than the threshold, then we have no match.
  if (neighbors.n_neighbors() != neighbors.max_neighbors()
      || neighbors.max_angular_distance() > max_angular_distance) {
    return false;
  }

  match = neighbors.nearest_neighbor();
  return true;
}

// Moving private method here since it's called by public methods following
void tree_union::get_node_level_pixels(const pixel& pix, pixel_vector* pixels) const {
  if (!pixels->empty())
    pixels->clear();

  // Assemble an array of pixels at the level of the tree_map or higher from
  // the input pixel.
  if (pix.level() >= level_) {
    pixels->push_back(pix);
  } else if (pix.level() < level_) {
    for (pixel c = pix.child_begin(level_); c != pix.child_end(level_); c
        = c.next()) {
      pixels->push_back(c);
    }
  }
}

long tree_union::n_points(const pixel& pix) const {
  if (!may_intersect(pix)) {
    return 0;
  }

  pixel_vector pixels;
  get_node_level_pixels(pix, &pixels);

  // Scan over those pixels and gather the aggregate to return.
  long total_points = 0;
  for (pixel_iterator iter = pixels.begin(); iter != pixels.end(); ++iter) {
    tree_map_iterator t_iter = tree_map_.find(iter->parent(level_).id());
    if (t_iter != tree_map_.end()) {
      total_points = t_iter->second->n_points(*iter);
    }
  }

  return total_points;
}

double tree_union::weight(const pixel& pix) const {
  if (!may_intersect(pix)) {
    return 0.0;
  }

  pixel_vector pixels;
  get_node_level_pixels(pix, &pixels);

  // Scan over those pixels and gather the aggregate to return.
  double total_weight = 0;
  for (pixel_iterator iter = pixels.begin(); iter != pixels.end(); ++iter) {
    tree_map_iterator t_iter = tree_map_.find(iter->parent(level_).id());
    if (t_iter != tree_map_.end()) {
      total_weight = t_iter->second->weight(*iter);
    }
  }

  return total_weight;
}

void tree_union::points(point_vector* points) const {
  if (!points->empty())
    points->clear();

  for (tree_map_iterator iter = tree_map_.begin(); iter != tree_map_.end(); ++iter) {
    point_vector tmp_points;
    iter->second->points(&tmp_points);

    for (point_iterator p_iter = tmp_points.begin(); p_iter
        != tmp_points.end(); ++p_iter) {
      points->push_back(*p_iter);
    }
  }
}

void tree_union::points(const pixel& pix, point_vector* points) const {
  if (!points->empty())
    points->clear();

  if (!may_intersect(pix)) {
    return;
  }

  pixel_vector pixels;
  get_node_level_pixels(pix, &pixels);

  // Scan over those pixels and gather the aggregate to return.
  for (pixel_iterator iter = pixels.begin(); iter != pixels.end(); ++iter) {
    tree_map_iterator t_iter = tree_map_.find(iter->parent(level_).id());
    if (t_iter != tree_map_.end()) {
      point_vector tmp_points;
      t_iter->second->points(*iter, &tmp_points);

      for (point_iterator p_iter = tmp_points.begin(); p_iter
          != tmp_points.end(); ++p_iter) {
        points->push_back(*p_iter);
      }
    }
  }
}

void tree_union::clear() {
  for (tree_map_iterator iter = tree_map_.begin(); iter != tree_map_.end(); ++iter) {
    iter->second->clear();
    delete iter->second;
  }
  tree_map_.clear();
  nodes_.clear();
  bound_.clear();
}

void tree_union::calculate_area() {
  if (!modified_) {
    return;
  }

  area_ = area();

  modified_ = false;
}

double tree_union::area() const {
  if (!modified_)
    return area_;

  double total_area = 0.0;
  for (tree_map_iterator iter = tree_map_.begin(); iter != tree_map_.end(); ++iter) {
    total_area += iter->second->covering_fraction()
        * iter->second->exact_area();
  }

  return total_area;
}

bool tree_union::contains(const point& p) const {
  return contains(p.to_pixel());
}

bool tree_union::contains(const pixel& pix) const {
  pixel_vector pixels;
  get_node_level_pixels(pix, &pixels);

  for (pixel_iterator iter = pixels.begin(); iter != pixels.end(); ++iter) {
    if (!may_intersect(*iter)) {
      return false;
    }
  }

  return true;
}

double tree_union::contained_area(const pixel& pix) const {
  if (!may_intersect(pix)) {
    return 0.0;
  }

  pixel_vector pixels;
  get_node_level_pixels(pix, &pixels);

  // Scan over those pixels and gather the aggregate to return.
  double total_area = 0.0;
  for (pixel_iterator iter = pixels.begin(); iter != pixels.end(); ++iter) {
    tree_map_iterator t_iter = tree_map_.find(iter->parent(level_).id());
    if (t_iter != tree_map_.end()) {
      total_area = iter->exact_area()
          * t_iter->second->covering_fraction(*iter);
    }
  }

  return total_area;
}

bool tree_union::may_intersect(const pixel& pix) const {
  pixel_set_iterator iter = nodes_.lower_bound(pix);
  if (iter != nodes_.end() && iter->range_min() <= pix.range_max())
    return true;
  return iter != nodes_.begin() && (--iter)->range_max() >= pix.range_min();
}

void tree_union::initialize_bound() {
  if (tree_map_.empty()) {
    bound_ = circle_bound();
  }

  initialized_bound_ = true;
  bound_ = get_bound();
}

point tree_union::get_center() const {
  if (initialized_bound_) {
    return bound_.axis();
  }

  double x = 0.0, y = 0.0, z = 0.0;
  for (tree_map_iterator iter = tree_map_.begin(); iter != tree_map_.end(); ++iter) {
    double area = iter->second->exact_area();
    point center = iter->second->center_point();
    x += area * center.unit_sphere_x();
    y += area * center.unit_sphere_y();
    z += area * center.unit_sphere_z();
  }

  return point(x, y, z, 1.0);
}

circle_bound tree_union::get_bound() const {
  if (initialized_bound_) {
    return bound_;
  }

  circle_bound bound = circle_bound(get_center(), 0.0);
  for (tree_map_iterator iter = tree_map_.begin(); iter != tree_map_.end(); ++iter) {
    bound.add_circle_bound(iter->second->get_bound());
  }

  return bound;
}

void tree_union::get_covering(pixel_vector* pixels) const {
  // this is the default mode for a pixel covering. For this as for other
  // coverings we use only 8 pixels at max to cover our union.
  get_covering(DEFAULT_COVERING_PIXELS, pixels);
}

void tree_union::get_covering(long max_pixels, pixel_vector* pixels) const {
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

void tree_union::get_covering(double fractional_area_tolerance,
    pixel_vector* pixels) const {
  if (!pixels->empty())
    pixels->clear();
  pixels->reserve(nodes_.size());

  for (tree_map_iterator iter = tree_map_.begin(); iter != tree_map_.end(); ++iter) {
    pixels->push_back(*(iter->second));
  }
}

void tree_union::get_interior_covering(int max_level, pixel_vector* pixels) const {
  if (!pixels->empty())
    pixels->clear();
  pixels->reserve(nodes_.size());

  for (tree_map_iterator iter = tree_map_.begin(); iter != tree_map_.end(); ++iter) {
    if (iter->second->level() <= max_level) {
      pixels->push_back(*(iter->second));
    }
  }
}

void tree_union::get_simple_covering(int level, pixel_vector* pixels) const {
  // For a simple covering of a pixel_union we don't need anything fancy like
  // FloodFill from S2 (at least for this method). If we want a simple covering
  // of a pixel_union then we just need to loop through the pixels, make parents
  // out of children
  if (!pixels->empty())
    pixels->clear();

  for (tree_map_iterator iter = tree_map_.begin(); iter != tree_map_.end(); ++iter) {
    tree_pixel* pix_iter = iter->second;
    if (pix_iter->level() < level) {
      for (pixel c = pix_iter->child_begin(); c != pix_iter->child_end(); c
          = c.next()) {
        if (!pixels->empty() && pixels->back().contains(c))
          continue;
        pixels->push_back(c);
      }
    } else if (pix_iter->level() == level) {
      if (!pixels->empty() && pixels->back().contains(*pix_iter))
        continue;
      pixels->push_back(*pix_iter);
    } else if (pix_iter->level() > level) {
      pixel pix = pix_iter->parent(level);
      if (!pixels->empty() && pixels->back().contains(pix))
        continue;
      pixels->push_back(pix);
    }
  }
}

void tree_union::get_center_covering(int level, pixel_vector* pixels) const {
  get_simple_covering(level, pixels);
}

} // end namespace s2omp
