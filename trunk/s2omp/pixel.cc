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
#include "angular_bin-inl.h"
#include "circle_bound.h"

bool pixel::edge_distances(point& p, double& near_edge_distance,
		double& far_edge_distance) {
	// First, find the nearest vertex.
	double near_edge_distance = 100.0;
	double far_edge_distance = -100.0;
	for (int k = 0; k < 4; k++) {
		double costheta = vertex(k).dot(p);
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
		edges.push_back(edge(k));
		edge_caps.push_back(circle_bound.from_height(edge(k), 0.0));
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
