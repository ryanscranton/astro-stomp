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
}

void tree_pixel::initialize_node(uint16_t max_points) {
	initialized_children_ = false;
	weight_ = 0.0;
	maximum_points_ = 0;
	point_count_ = 0;
	maximum_points_ = max_points;
	point_count_ = 0;
}

explicit tree_pixel::tree_pixel(uint64 id) {
	set_id(id);
	initialize_node(kDefaultMaxPoints);
}

tree_pixel::tree_pixel(uint64 id, uint16_t max_points) {
	set_id(id);
	initialize_node(max_points);
}

tree_pixel::~treel_pixel() {
	clear();
}

tree_pixel* tree_pixel::from_point(const point& p, int level,
		uint16_t max_points) {
	return new tree_pixel(point_to_id(p, level), max_points);
}

tree_pixel* tree_pixel::from_pixel(const pixel& pix, uint16_t max_points) {
	return new tree_pixel(pix.id(), max_points);
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

double tree_pixel::find_weighted_pairs(const annulus_bound& bound) const {
	double pair_weight = 0;

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
					pair_weight += (*iter)->find_weighted_pairs(bound);
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
			covered_area += total_area;
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
		// case 2, this tree pixel contains pix and has points
		if (!initialized_children_) {
			// case 2.1, pixel is contained but tree_pixel has no children
			// we then need to loop over all the points and see if pixel contains any
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

void tree_pixel::points(point_vector* point_vec) const {
	if (!point_vec->empty()) point_vec->clear();

	if (point_count_ > 0) {
		point_vec->reserve(point_count_);
		if (!initialized_children_) {
			for (point_ptr_iterator iter = points_.begin();
					iter != points_.end(); iter++) {
				point_vec->push_back(*(*iter));
			}
		} else if (initialized_children_) {
			for (tree_ptr_iterator iter = children_.begin();
					iter != children_.end(); iter++) {
				point_vector tmp_points;
				(*iter)->points(&tmp_points);
				for (point_iterator p_iter = tmp_points.begin();
						p_iter != tmp_points.end(); p_iter++) {
					point_vec->push_back(*p_iter);
				}
			}
		}
	}
}

void tree_pixel::points(pixel& pix, point_vector* point_vec) const {
	if (!point_vec->empty()) point_vec->clear();

	if (point_count_ > 0) {
		point_vec->reserve(point_count_);
		if (pix.contains(*this)) {
			points(point_vec);
		} else if (contains(pix)) {
			if (!initialized_children_) {
				for (point_ptr_iterator iter = points_.begin();
						iter != points_.end(); iter++) {
					if (pix.contains(*(*iter))) point_vec->push_back(*(*iter));
				}
			} else if (initialized_children_) {
				for (tree_ptr_iterator iter = children_.begin();
						iter != children_.end(); iter++) {
					point_vector tmp_points;
					(*iter)->points(&tmp_points);
					for (point_iterator p_iter = tmp_points.begin();
							p_iter != tmp_points.end(); p_iter++) {
						point_vec.push_back(*p_iter);
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
				iter != points_.end(); iter++) {
			delete *iter;
		}
	}
	points_.clear();

	if (initialized_children_) {
		for (tree_ptr_iterator iter = children_.begin();
				iter != children_.end(); ++iter) {
			iter->clear();
			delete *iter;
		}
	}
	children_.clear();
}

bool tree_pixel::initialize_children() {
	initialized_children_ = false;
	if (!is_leaf()) {
		children_.reserve(4);

		pixel_vector* tmp_children;
		children(tmp_children);
		for (pixel_iterator iter = tmp_children->begin();
				iter != tmp_children->end(); ++iter) {
			tree_pixel t_pix =
					new tree_pixel(iter->id(), maximum_points_);
			children_.push_back(t_pix);
		}
		initialized_children_ = true;
		/*
		bool transferred_point_to_children = false;
		for (point_ptr_iterator p_iter = points_.begin();
				p_iter != points_.end(); ++p_iter) {
			transferred_point_to_children = false;
			for (tree_ptr_iterator iter = children_.begin();
					iter != children_.end(); ++iter) {
				if ((*iter)->add_point(p_iter))
					// pretty sure I copied this over the same but it looks like there is
					// a bug here. If an early point fails and a later point succeeds
					// then the output of this loop will still be true. Look below for
					// better version?
					transferred_point_to_children = true;
			}
		}*/
		bool transferred_points_to_children = true;
		point_ptr_iterator p_iter = points_.begin();
		while (p_iter != points_.end() && transferred_points_to_children) {
			transferred_points_to_children = false;
			for (tree_ptr_iterator iter = children_.begin();
					iter != children_.end(); ++iter) {
				if ((*iter)->add_point(p_iter))
					transferred_points_to_children = true;
			}
			++p_iter;
		}
		if (!transferred_points_to_children) initialized_children_ = false;
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

void tree_pixel::neighbor_recursion(point& p, tree_neighbor& neighbor) {
	neighbor.add_node();

	if (!initialized_children_) {
		// We have no sub-nodes in this tree, so we'll just iterate over the
		// points here and take the nearest N neighbors.
		for (point_ptr_iterator iter = points_.begin();
				iter != points_.end(); ++iter) {
			neighbor.test_point(*iter);
		}
	} else {
		// This node is the root node for our tree, so we first find the sub-node
		// that contains the point and start recursing there.
		//
		// While we iterate through the nodes, we'll also calculate the edge
		// distances for those nodes that don't contain the point and store them
		// in a priority queue.  This will let us do a follow-up check on nodes in
		// the most productive order.
		pixel_queue pix_queue;
		for (tree_ptr_iterator iter = children_.begin();
				iter != children_.end(); ++iter) {
			if ((*iter)->contains(p)) {
				(*iter)->neighbor_recursion(ang, neighbor);
			} else {
				distance_pixel_pair dist_pair((*iter)->nearest_edge_distance(p),
						(*iter));
				pix_queue.push(dist_pair);
			}
		}

		// That should give us back a TreeNeighbor object that contains a workable
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
			if (pix_distance < neighbor.max_distance()) {
				pix_iter->neighbor_recursion(p, neighbor);
			}
			pix_queue.pop();
		}
	}
}

tree_neighbor::tree_neighbor(const point& reference_point) {
	tree_neighbor(point, 1, 100.0);
}

tree_neighbor::tree_neighbor(const point& reference_point,
		uint8_t n_neighbors) {
	tree_neighbor(reference_point_, n_neighbors, 100.0);
}

tree_neighbor::tree_neighbor(const point& reference_point, uint8_t n_neighbors,
			double max_distance) {
	reference_point_ = reference_point;
	n_neighbors_ = n_neighbors;
	max_distance_ = max_distance;
	n_nodes_visited_ = 0;
}

void tree_neighbor::nearest_neighbors(const point_vector& p_vect,
		bool save_neighbors) {
	if (!p_vect.empty()) p_vect.clear();
	std::vector<distance_point_pair> backup_copy;

	while(!point_queue_.empty()) {
		distance_point_pair dist_pair = point_queue_.top();
		point_queue_.pop();

		point tmp_point(dist_pair.second->unit_sphere_x(),
				dist_pair.second->unit_sphere_y(),
				dist_pair.second->unit_sphere_z(),
				dist_pair.second->weight());

		p_vect.push_back(tmp_point);
		backup_copy.push_back(dist_pair);
	}

	if (save_neighbors) {
		for (uint8_t i=0; i < backup_copy.size(); i++) {
			point_queue_.push(backup_copy[i]);
		}
	}
}

bool tree_neighbor::test_point(point* test_point) {
	bool kept_point = false;
	double costheta = reference_point_.dot(test_point);
	double sin2theta = 1.0 - costheta*costheta;

	if (sin2theta < max_distance_ || n_neighbors() < max_neighbors()) {
		kept_point = true;
		if (n_neighbors() == max_neighbors) point_queue_.pop();
		distance_point_pair dist_pair(sin2theta, test_point);
		point_queue_.push(dist_pair);
		max_distance_ = point_queue_.top().first;
	}
	return kept_point;
}

double tree_neighbor::max_angular_distance() {
	return RAD_TO_DEG*asin(sqrt(fabs(max_distance_)));
}

} //end namespace s2omp

