/*
 * pixel.cc
 *
 *  Created on: Jul 10, 2012
 *      Author: cbmorrison
 */




double pixel::nearest_edge_distance(point& p) {
	point edge0 = edge(0);
	point edge1 = edge(1);
	point edge2 = edge(2);
	point edge3 = edge(3);

	// test if point is within any edges using cap_bounds. If it is compute the
	// closest approach to the other edges.

	// one is not the first edge just the first one we test. I'm not really
	// coding this right now. Just putting down the equation.

	// the basic idea of this is that we need to find normal for the great circle
	// that is such that it goes through the point p and is normal to that of the
	// edge. The easiest way to do this is think in terms of intersecting planes.
	// if the edge_1 is the normal defining the edge great circle then we need to
	// solve the equations.

	// edge_1.dot(p.cross(x)) == 0 and edge_1.dot(x) == 0

	// where x is the point we are solving for. The last constraint is that
	// x.dot(x) == 1. The solution in terms of cross products is shown below.
	point cross_prod = edge_1.cross(p);

	double norm = edge_1.unit_sphere_y()*cross_prod.unit_sphere_x() -
			edge_1.unit_sphere_x()*cross_prod.unit_sphere_y();
	double x = (edge_1.unit_sphere_y()*cross_prod.unit_sphere_z() +
			edge_1.unit_sphere_z()*cross_prod.unit_sphere_x())/norm;
	double y = (edge_1.unit_sphere_x()*cross_prod.unit_sphere_z() +
				edge_1.unit_sphere_z()*cross_prod.unit_sphere_x())/norm;
	double z = 1.0;
	norm = sqrt(x*x + y*y + z*z);
	point closest_canidate = point(x/norm, y/norm, z/norm, 1.0);
	// now have one of the closest approach canidates. The other looks the same.

	return distance;
}
