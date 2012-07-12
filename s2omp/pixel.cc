/*
 * pixel.cc
 *
 *  Created on: Jul 10, 2012
 *      Author: cbmorrison
 */




double pixel::nearest_edge_distance(point& p) {
	// we need the edges for calculating if the point is between two of the edges
	// as well as defining the circle bounds
	point edge0 = edge(0);
	point edge1 = edge(1);
	point edge2 = edge(2);
	point edge3 = edge(3);

	// from these edges we create hemisphere cirlce_bounds
	circle_bound cap0 = circle_bound.from_height(edge0, 0.0);
	circle_bound cap1 = circle_bound.from_height(edge1, 0.0);
	circle_bound cap2 = circle_bound.from_height(edge2, 0.0);
	circle_bound cap3 = circle_bound.from_height(edge3, 0.0);

	// We're going to possibly be testing contains on the bounds above several
	// times. As such it would be useful to store the boolean after it is
	// determined.
	bool in_cap0, in_cap2 = false;

	// we have to do these tests for any of the conditions below, so we run them
	// first and save the results.
	in_cap0 = cap0.contains(p);
	in_cap2 = cap2.contains(p);

	// test if the point is contained between edge 0 and 2. If it is the nearest
	// point is on either edge 1 or 3.
	if (in_cap0 && in_cap2) {
		// test if point is inside edge3's cap. If it is then the nearest point
		// is on edge1.
		if (cap3.contains(p)) {
			// test if point is within any edges using cap_bounds. If it is compute
			// the closest approach to the other edges.

			// the basic idea of this is that we need to find the normal for the great
			// circle that goes through the point p and is normal to that
			// of the edge. The easiest way to do this is think in terms of
			// intersecting planes. if the nearest_edge is the normal defining the
			// edge great circle then we need to solve the equations.

			// edge_1.dot(p.cross(x)) == 0 and edge_1.dot(x) == 0

			// where x is the point we are solving for. The last constraint is that
			// x.dot(x) == 1. The solution in terms of cross products is shown below.
			nearest_point = edge1.cross(p.cross(edge1)); // need to normalize
			distance = p.angular_distance(nearest_point);
		} else if (cap1.contains(p)) {
			// same as above, different index.
			nearest_point = edge3.cross(p.cross(edge3));
			distance = p.angular_distance(p);
		}
	} else if (cap1.contains(p) && cap3.contains(p)) {
		// same as above, different indexes.
		if (in_cap2) {
			nearest_point = edge0.cross(p.cross(edge0)); // need to normalize
			distance = p.angular_distance(nearest_point);
		} else if (in_cap0) {
			nearest_point = edge2.cross(p.cross(edge2)); // need to normalize
			distance = p.angular_distance(nearest_point);
		}
	} else {
		// the nearest distance has to be one of the vertexes check the bounds
		// this time if the point is contained by both bounds, the closest approach
		// is the vertex on the other side.
		if (in_cap0 && cap1.contains(p)) {
			distance = p.angular_distance(vertex(1));
		} else if (cap1.contains(p) && in_cap2) {
			distance = p.angular_distance(vertex(2));
		} else if (in_cap2 && cap3.contains(p)) {
			distance = p.angular_distance(vertex(3));
		} else {
			distance = p.angular_distance(vertex(0));
		}
	}

	return distance;
}
