/*
 * tree_union.cc
 *
 *  Created on: Jul 16, 2012
 *      Author: cbmorrison
 */

#include "tree_union.h"

explicit tree_union::tree_union(int level) {
  tree_union(level, 200);
}

tree_union::tree_union(int level, max_points) {
  maximum_points_ = max_points;
  nodes_ = 0;
  point_count_ = 0;
  level_ = level;
  weight_ = 0.0;
  area_ = 0.0;
  modified_ = false;
}

bool tree_union::add_point(point& p) {
  tree_map_iterator iter = tree_map_.find(p.id(level_));

  if (iter == tree_map_.end()) {
    // This point initializes a new pixel which we need to initialize.
    tree_map_.insert(std::pair<uint64, tree_pixel*>(p.id(level_),
        tree_pixel.from_point(p, level_, maximum_points_)));
    iter = tree_map_.find(p.id());
    if (iter == tree_map_.end()) {
      std::cout << "s2omp::tree_union::add_point - " <<
          "Creating new tree_union node failed. Exiting.\n";
      exit(2);
    }
    area += (*iter).second->exact_area();
  }
  bool added_point = (*iter).second->add_point(p);
  if (added_point) {
    point_count_++;
    weight += p.weight();
  }
}

uint32_t tree_union::find_pairs(const annulus_bound& bound) const {
  uint32_t pairs = 0;
  pixel_vector* pixels;
  bound.simple_covering(level_, pixels);

  for (pixel_iterator iter = pixels->begin(); iter != pixels->end(); ++iter) {
    tree_map_iterator tree_iter = tree_map_.find(iter->id());
    if (tree_iter == tree_map_.end()) continue;
    pairs += (*tree_iter).second->find_pairs(bound);
  }
  pixels->clear();
  delete pixels;
  return pairs;
}

void tree_union::find_pairs(const point_vector& p, angular_bin* bin) const {
  // Loop through points and find the number of pairs within the annulus defined
  // by angular_bin. Store values in bin.
}

double tree_union::find_weighted_pairs(const annulus_bound& bound) const {
  double wpairs = 0.0;
  pixel_vector* pixels;
  bound.simple_covering(level_, pixels);
  for (pixel_iterator iter = pixels->begin(); iter != pixels->end(); ++iter) {
    tree_map_iterator tree_iter = tree_map_.find(iter->id());
    if (tree_iter == tree_map_.end()) continue;
    wpairs += (*tree_iter).second->find_weight_pairs(bound);
  }
  pixel->clear();
  delete pixels;
  return wpairs;
}

void tree_union::find_pairs(const point_vector& p, angular_bin* bin) const {
  // Loop through points and find the number of weighted pairs within the
  // annulus defined by angular_bin. Store values in bin.
}

