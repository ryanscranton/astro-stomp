/*
 * circle_bound.h
 *
 *  Created on: Jul 5, 2012
 *      Author: scranton
 */

#ifndef CIRCLE_BOUND_H_
#define CIRCLE_BOUND_H_

#include <stdint.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "core.h"
#include "geometric_bound_interface.h"
#include "MersenneTwister.h"

namespace s2omp {
class point;  // class declaration in stomp_angular_coordinate.h
class angular_bin;         // class declaration in stomp_angular_bin.h
class circle_bound;

typedef std::vector<circle_bound> circle_vector;
typedef circle_vector::iterator circle_iterator;

class circle_bound: public geometric_bound {
public:
	virtual ~circle_bound();

	static circle_bound* from_angular_bin(const point& axis,
			const angular_bin bin);
	static circle_bound* from_radius(const point& axis, double radius_degrees);
	static circle_bound* from_height(const point& axis, double height);

	inline point axis() {return axis_;}
	inline double radius() {return radius_;}
	inline double height() {return height_;}

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
	circle_bound();
	void initialize_random();

	S2::S2Cap cap_;
	point axis_;
  bool initialized_random_;
	double radius_, height_;

	//These variables below are for generating random points and will not be
	// initialized unti initialize_random() is called.
	point great_cirlce_norm_;
	double rotate_;
	bool initialized_random_;
	MTRand mtrand;
};

} // end namespace s2omp


#endif /* CIRCLE_BOUND_H_ */
