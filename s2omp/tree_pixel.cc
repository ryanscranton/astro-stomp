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

bool tree_pixel::add_point(const point* p) {

	bool added_to_pixel = false;
	if (contains(*p)) {
		if (point_count_ < maximum_points_ || is_leaf()) {
			if (point_count_ == 0)
				points_.reserve(maximum_points_);
			points_.push_back(p);
			added_to_pixel = true;
		} else {
			if (!initialized_children_) {
				if (!initialize_children()) {
					std::cout << "s2omp::tree_pixel::add_point - "
							<< "Failed to initialize children.  Exiting.\n";
					exit(2);
				}
			} // end initialize if
			for (tree_ptr_iterator iter = children_.begin();
					iter != children_.end(); iter++) {
				if ((*iter)->contains(*p)) {
					added_to_pixel = (*iter)->add_point(p);
					if (added_to_pixel) points_.push_back(p);
					break;
				}
			} //end children loop
		} //end if(point_count_)
	} else {
		added_to_pixel = false;
	} //end contains if

	if (added_to_pixel) {
		weight_ += p->weight();
		point_count_++;
	}

	return added_to_pixel;
}

bool tree_pixel::add_point(const point& p) {
	point* point_copy =
	    new point(point.unit_sphere_x(), point.unit_sphere_y(),
					  point.unit_sphere_z(), point.weight());
	return add_point(point_copy);
}

uint32_t tree_pixel::find_pairs(const annulus_bound& bound) const {
	uint32_t pair_count = 0;

	if (bound.may_intersect(*this)) {
		if (bound.contains(*this)) {
			pair_count += point_count_;
		} else {
			if (!initialized_children_) {
				for (point_ptr_iterator iter = points_.begin();
						iter != points_.end(); iter++) {
					if (bound.contains(*(*iter))) pair_count++;
				}
			} else {
				for (tree_ptr_iterator iter = children_.begin();
						iter != children_.end(); iter++) {
					pair_count += (*iter)->find_pairs(bound);
				}
			}
		}
	} else {
		// no intersection
		pair_count = 0;
	}
	return pair_count;
}

uint32_t tree_pixel::find_pairs(const annulus_bound& bound) const {
	uint32_t pair_weight = 0;

	if (bound.may_intersect(*this)) {
		if (bound.contains(*this)) {
			pair_weight += point_count_;
		} else {
			if (!initialized_children_) {
				for (point_ptr_iterator iter = points_.begin();
						iter != points_.end(); iter++) {
					if (bound.contains(*(*iter))) pair_weight += (*iter)->weight();
				}
			} else {
				for (tree_ptr_iterator iter = children_.begin();
						iter != children_.end(); iter++) {
					pair_weight += (*iter)->find_pairs(bound);
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

uint16_t tree_pixel::find_nearest_neighbor(const point& p,
		point* neighbor) const {
	point_vector p_vector;

	uint16_t nodes_visited = find_k_nearest_neighbors(p, 1, p_vector);

	neighbor = p_vector[0];

	return nodes_visited;
}

double tree_pixel::k_nearest_neighbor_distance(const point& p,
		uint8_t n_neighbors, uint16_t& nodes_visited) const {

	tree_neighbor neighbors(p, n_neighbors);

	neighbor_recursion(p, neighbors);

	nodes_visited = neighbors.nodes_visited();

	return neighbors.max_angular_distance();
}

double tree_pixel::nearest_neighbor_distance(const point& p,
		uint16_t& nodes_visited) const {

	return k_nearest_neighbor_distance(p, 1, nodes_visited);
}

bool tree_pixel::closest_match(const point& p, double max_angular_distance,
		point* match) const {

	tree_neighbor neighbors(p, 1, max_angular_distance);

	neighbor_recursion(p, neighbors);

	bool found_match = false;
	if (neighbors.neighbors() == neighbors.max_neighbors() &&
			neighbors.max_angular_distance() < max_angular_distance) {
		found_match = true;

		point_vector p_vector;
		neighbors.nearest_neighbors(p_vector, false);
		match = p_vector[0];
	}

	return found_match;
}

uint32_t tree_pixel::n_points() const {
	return point_count_;
}

double tree_pixel::weight() const {
	return weight_;
}

uint32_t tree_pixel::n_points(const pixel& pix) const {
	uint32_t contained_points = 0;

	if (pix.contains(*this)) {
		contained_points += point_count_;
	} else if (contains(pix)) {
		if (!initialized_children_) {
			for (point_ptr_iterator iter = points_.begin();
					iter != points_.end(); iter++) {
				if (pix.contains(*(*iter))) contained_points++;
			}
		} else if (initialized_children_) {
			for (tree_ptr_iterator iter = children_.begin();
					iter != children_.end(); iter++) {
				contained_points += (*iter)->n_points(pix);
			}
		}
	}
	return contained_points;
}

uint32_t tree_pixel::weight(const pixel& pix) const {
	uint32_t contained_weight = 0;

	if (pix.contains(*this)) {
		contained_weight += weight_;
	} else if (contains(pix)) {
		if (!initialized_children_) {
			for (point_ptr_iterator iter = points_.begin();
					iter != points_.end(); iter++) {
				if (pix.contains(*(*iter))) contained_weight += (*iter)->weight();
			}
		} else if (initialized_children_) {
			for (tree_ptr_iterator iter = children_.begin();
					iter != children_.end(); iter++) {
				contained_weight += (*iter)->weight(pix);
			}
		}
	}
	return contained_weight;
}

double tree_pixel::covering_fraction() const {
	double total_area = exact_area();
	double covered_area = 0.0;

	if (point_count_ > 0) {
		if (!initialized_children_) {
			covered_area += exact_area();
		} else if (initialized_children_) {
			for (tree_ptr_iterator iter = children_.begin();
					iter != children_.end(); iter++) {
				covered_area += (*iter)->covering_fraction()*
						(*iter)->exact_area();
			}
		}
	}
	return covered_area/total_area;
}

double tree_pixel::covering_fraction(pixel& pix) const {
	double covered_fraction = 0.0;

	if (pix.contains(*this)) {
		// case 1, pix is larger than the current tree pixel (or equal to)
		covered_fraction += exact_area()*covering_fraction()/pix.exact_area();
	} else if (contains(pix) && point_count_ > 0) {
		//case 2, this tree pixel contains pix and has points
		if (!initialized_children_) {
			//case 2.1, pixel is contained but tree_pixel has no children
			//we then need to loop over all the points and see if pixel contains any
			for (point_ptr_iterator iter = points_.begin();
					iter != points_.end(); iter++) {
				if (pix.contains(*(*iter))) {
					covered_fraction = 1.0;
					break;
				}
			}
		} else {// case 2.2 tree_pixel has children, need to loop over them
			for (tree_ptr_iterator iter = children_.begin();
					iter != children_.end(); iter++) {
				covered_fraction += (*iter)->covering_fraction(pix);
			}
		}
	}
	return covered_fraction;
}

void tree_pixel::points(point_vector& p_vect) const {
	if (!p_vect.empty()) p_vect.clear();
	p_vect.reserve(point_count_);

	if (!initialized_children_) {
		for (point_ptr_iterator iter = points_.begin();
				iter != points_.end(); iter++) {
			p_vect.push_back(*(*iter));
		}
	} else {
		for (tree_ptr_iterator iter = children_.begin();
				iter != children_.end(); iter++) {
			point_vector tmp_p_vect;
			(*iter)->points(tmp_p_vect);
			for (point_iterator p_iter = tmp_p_vect.begin();
					p_iter != tmp_p_vect.end(); p_iter++) {
				p_vect.push_back(*p_iter);
			}
		}
	}
}

void tree_pixel::points(pixel& pix, point_vector& p_vect) const {
	if (!p_vect.empty()) p_vect.clear();
	p_vect.reserve(point_count_);

	if (point_count_ > 0) {
		if (pix.contains(*this)) {
			points(p_vect);
		} else if (contains(pix)) {
			if (!initialized_children_) {
				for (point_ptr_iterator iter = points_.begin();
						iter != points_.end(); iter++) {
					if (pix.contains(*(*iter))) p_vect.push_back(*(*iter));
				}
			} else if (initialized_children_) {
				for (tree_ptr_iterator iter = children_.begin();
						iter != children_.end(); iter++) {
					point_vector tmp_p_vect;
					(*iter)->points(tmp_p_vect);
					for (point_iterator p_iter = tmp_p_vect.begin();
							p_iter != tmp_p_vect.end(); p_iter++) {
						p_vect.push_back(*p_iter);
					}
				}
			}
		}
	}
}

uint16_t tree_pixel::n_nodes() const {
	uint16_t total_nodes = 1;

	if (initialized_children_) {
		for (tree_ptr_iterator iter = children_.begin();
				iter != children_.end(); iter++) {
			total_nodes += (*iter)->n_nodes();
		}
	}
	return total_nodes;
}

void tree_pixel::clear() {
	if (!points_.empty()) {
		for (point_ptr_iterator iter = points_.begin();
				iter != points_.end(); iter++)
			delete *iter;
	}
	points_.clear();
	if (initialized_children_) {
		for (tree_ptr_iterator iter = children_.begin();
				iter != children_.end(); ++iter)
			delete *iter;
	}
	children_.clear();
}

bool tree_pixel::initialize_children() {
	initialized_children_ = false;
	if (!is_leaf()) {
		pixel_vector* tmp_children;
		children_.reserve(4);
		children(tmp_children);

		for (pixel_iterator iter = tmp_children->begin();
				iter != tmp_children->end(); ++iter) {
			tree_pixel t_pix =
					new tree_pixel(iter->id(), maximum_points_);
			children_.push_back(t_pix);
		}
		initialized_children_ = true;

		bool transferred_point_to_children = false;
		for (point_ptr_iterator p_iter = points_.begin();
				p_iter != points_.end(); ++p_iter) {
			transferred_point_to_children = false;
			for (tree_ptr_iterator iter = children_.begin();
					iter != children_.end(); ++iter) {
				if ((*iter)->add_point(p_iter))
					transferred_point_to_children = true;
			}
		}
		if (!transferred_point_to_children) initialized_children_ = false;
	} // end is_leaf check
	return initialized_children_;
}

uint32_t tree_pixel::direct_pair_count(annulus_bound& bound) {
	uint32_t pair_count = 0;

	for (point_ptr_iterator iter = points_.begin();
			iter != points_.end(); ++iter) {
		if (bound.contains(*(*iter))) pair_count++;
	}
	return pair_count;
}

double tree_pixel::direct_weighted_pairs(annulus_bound& bound) {
	double pair_weight = 0;

	for (point_ptr_iterator iter = points_.begin();
			iter != points_.end(); ++iter) {
		if (bound.contains(*(*iter))) pair_weight += (*iter)->weight;
	}
	return pair_weight;
}

/* Don't have all the methods neede for this.
void tree_pixel::neighbor_recursion(point& p, tree_neighbor& neighbor) {
	neighbor.add_node();

	if (!points_.empty()) {
		for (point_ptr_iterator iter = points_.begin();
				iter != points_.end(); ++iter) {
			neighbor.test_point(*iter);
		}
	} else {
		pixel_queue pix_queue;
		for (tree_ptr_iterator iter = children_.begin();
				iter != children_.end(); ++iter) {
			if ((*iter)->contains(p)) {
				(*iter)->neighbor_recursion(ang, neighbor);
			} else {
				double min_edge_distance, max_edge_distance;
				(*iter)->
			}
		}
	}
}
*/

} //end namespace s2omp

