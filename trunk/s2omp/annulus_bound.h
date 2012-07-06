/*
 * annulus_bound.h
 *
 *  Created on: Jul 5, 2012
 *      Author: scranton
 */

#ifndef ANNULUS_BOUND_H_
#define ANNULUS_BOUND_H_

namespace s2omp {

class annulus_bound: public geometric_bound {
public:
	virtual ~annlus_bound();

	static annulus_bound* from_radii(const point& axis,
			double inner_theta_degrees, double outer_theta_degrees);
	static annulus_bound* from_cosradii(const point& axis,
			double inner_costheta, double outer_costheta);

	// API from geometric_bound
	virtual bool is_empty();
	virtual long size();
	virtual void clear();
	virtual void area();

	virtual bool contains(const point& p);
	virtual bool contains(const pixel& pix);
	virtual bool contains(const geometric_bound& b);

	virtual double contained_area(const pixel& pix);
	virtual bool may_intersect(const pixel& pix);

	virtual void covering(pixel_vector* pixels);
	virtual void covering(int max_pixels, pixel_vector* pixels);
	virtual void covering(double fractional_area_tolerance, pixel_vector* pixels);
	virtual void simple_covering(int level, pixel_vector* pixels);

	virtual circle_bound get_bound();

	virtual point get_random_point();
	virtual void get_random_points(long n_points, pixel_vector* points);

private:
	annulus_bound();
	S2Cap inner_cap_, outer_cap_;
	point axis_;
	double inner_radius_, outer_radius_;
};

} // end namespace s2omp

#endif /* ANNULUS_BOUND_H_ */
