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
	initialized_children_ = false;
	weight_ = 0.0;
	maximum_points_ = 0;
	point_count_ = 0;
}

explicit tree_pixel::tree_pixel(uint64 id) {
	tree_pixel(id, 200);
}

tree_pixel::tree_pixel(uint64 id, uint16_t max_points) {

	set_id(id);
	initialized_children_ = false;
	weight_ = 0.0;
	maximum_points_ = 0;
	point_count_ = 0;
	maximum_points_ = max_points;
	point_count_ = 0;
}

tree_pixel tree_pixel::from_point(const point& p, int level,
		uint16_t max_points) {

	uint64 id = point_to_id(p, level);
	tree_pixel t_pix(id, max_points);
	return t_pix;
}

tree _pixel tree_pixel::from_pixel(const pixel& pix, uint16_t max_points) {

	tree_pixel t_pix(pix.id(), max_points);
	return t_pix;
}

bool tree_pixel::add_point(const point& p) {

	bool added_to_pixel = false;
	if (contains(p)) {
		if (point_count_ < maximum_points_ || is_leaf()) {
			if (point_count_ == 0)
				points_->reserve(maximum_points_);
			points_->push_back(p);
			added_to_pixel = true;
		} else {
			if (!initialized_children_) {
				if (!initialize_children()) {
					std::cout << "s2omp::tree_pixel::add_point - "
							<< "Failed to initialize children.  Exiting.\n";
					exit(2);
				}
			} // end initialize if
			for (tree_ptr_iterator t_ptr_iter = children_->begin();
					t_ptr_iter != children_->end(); t_ptr_iter++) {
				if (t_ptr_iter->contains(p)) {
					added_to_pixel = t_ptr_iter->add_point(p);
					break;
				}
			} //end children loop
		} //end if(point_count_)
	} else {
		added_to_pixel = false;
	} //end contains if

	if (added_to_pixel) {
		add_to_weight(p.weight());
		point_count_++;
	}

	return added_to_pixel;
}

uint32_t tree_pixel::find_pairs(const annulus_bound& bound) const {
	uint32_t pair_count = 0;

	if (bound.may_intersect(this)) {
		if (bound.contains(this)) {
			pair_count += point_count_;
		} else {
			if (!initialized_children_ || is_leaf()) {
				for (point_ptr_iterator p_iter = points_->begin();
						p_iter != points_->end(); p_iter++) {
					if (bound.contains(*p_iter)) pair_count++;
				}
			} else {
				for (tree_ptr_iterator t_ptr_iter = children_->begin();
						t_ptr_iter != children_->end(); t_ptr_iter++) {
					pair_count += t_ptr_iter->find_pairs(bound);
				}
			}
		}
	} else {
		// no intersection
		pair_count = 0;
	}
	return pair_cout;
}

uint32_t tree_pixel::find_pairs(const annulus_bound& bound) const {
	uint32_t pair_weight = 0;

	if (bound.may_intersect(this)) {
		if (bound.contains(this)) {
			pair_weight += point_count_;
		} else {
			if (!initialized_children_ || is_leaf()) {
				for (point_ptr_iterator p_iter = points_->begin();
						p_iter != points_->end(); p_iter++) {
					if (bound.contains(*p_iter)) pair_weight += p_iter->weight();
				}
			} else {
				for (tree_ptr_iterator t_ptr_iter = children_->begin();
						t_ptr_iter != children_->end(); t_ptr_iter++) {
					pair_weight += t_ptr_iter->find_pairs(bound);
				}
			}
		}
	} else {
		// no intersection
		pair_weight = 0;
	}
	return pair_weight;
}

uint16_t find_k_nearest_neighbors(const point& p, uint8_t n_neighbors,
			point_vector* neighbor_vector) const {
	tree_neighbor neighbors(p, n_neighbors);

	neighbor_recursion(p, neighbors);

	neighbors.nearest_neighbors(neighbor_vector, false);

	return neighbors.nodes_visited();
}

uint32_t tree_pixel::n_points() const {
	return point_count_;
}

double tree_pixel::weight() const {
	return weight_;
}

uint32_t tree_pixel::n_points(const pixel& pix) const {
	uint32_t contained_points = 0;

	if (pix.contains(this)) {
		contained_points += point_count_;
	} else if (contains(pix)) {
		if (!initialized_children_ || is_leaf()) {
			for (point_ptr_iterator p_iter = points_->begin();
					p_iter != points_->end(); p_iter++) {
				if (pix.contains(*p_iter)) contained_points++;
			}
		} else if (initialized_children_) {
			for (tree_ptr_iterator t_ptr_iter = children_->begin();
					t_ptr_iter != children_->end(); t_ptr_iter++) {
				contained_points += t_ptr_iter->n_points(pix);
			}
		}
	}
	return contained_points;
}

uint32_t tree_pixel::weight(const pixel& pix) const {
	uint32_t contained_weight = 0;

	if (pix.contains(this)) {
		contained_weight += weight_;
	} else if (contains(pix)) {
		if (!initialized_children_ || is_leaf()) {
			for (point_ptr_iterator p_iter = points_->begin();
					p_iter != points_->end(); p_iter++) {
				if (pix.contains(*p_iter)) contained_weight += p_iter->weight();
			}
		} else if (initialized_children_) {
			for (tree_ptr_iterator t_ptr_iter = children_->begin();
					t_ptr_iter != children_->end(); t_ptr_iter++) {
				contained_weight += t_ptr_iter->n_points(pix);
			}
		}
	}
	return contained_weight;
}

bool tree_pixel::initialize_children() {
	initialized_children_ = false;
	if (!is_leaf()) {

	} // end is_leaf check
}

} //end namespace s2omp

