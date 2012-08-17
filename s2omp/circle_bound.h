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
#include "pixel.h"
#include "point.h"
#include "bound_interface.h"
#include "MersenneTwister.h"

namespace s2omp {

typedef std::vector<circle_bound> circle_vector;
typedef circle_vector::iterator circle_iterator;

class pixel;
class point;  // class declaration in stomp_angular_coordinate.h
class angular_bin;         // class declaration in stomp_angular_bin.h
class bound_interface;

class circle_bound;

class circle_bound: public bound_interface {
public:
	virtual ~circle_bound();

	static circle_bound* from_angular_bin(const point& axis,
			const angular_bin& bin);
	static circle_bound* from_radius(const point& axis, double radius_degrees);
	static circle_bound* from_height(const point& axis, double height);

	inline point axis() {return axis_;}
	inline double radius() {return acos(1 - height_)*DEG_TO_RAD;}
	inline double height() {return height_;}

	// API from geometric_bound
	virtual bool is_empty();
	virtual long size();
	virtual void clear();
	virtual double area();

	virtual bool contains(const point& p);
	virtual bool contains(const pixel& pix);
	virtual double contained_area(const pixel& pix);

	virtual bool may_intersect(const pixel& pix);

	virtual void covering(pixel_vector* pixels);
	virtual void covering(int max_pixels, pixel_vector* pixels);
	virtual void covering(double fractional_area_tolerance, pixel_vector* pixels);
	virtual void interior_covering(int max_level, pixel_vector* pixels);
	virtual void simple_covering(int level, pixel_vector* pixels);

	virtual circle_bound get_bound();

	virtual point get_random_point();
	virtual void get_random_points(long n_points, point_vector* points);
	virtual point get_weighted_random_point(const point_vector& points);
	virtual void get_weighted_random_points(long n_points, point_vector* points,
	    const point_vector& input_points);

private:
	circle_bound();
	circle_bound(const point& axis, double height);

	bool intersects(const pixel& pix, const point_vector& vertices);
	circle_bound* complement();

	point axis_;
	double height_;

	//These variables below are for generating random points and will not be
	// initialized until initialize_random() is called.
	void initialize_random();
	point great_circle_norm_;
	double rotate_;
	bool initialized_random_;
	MTRand mtrand;
};

} // end namespace s2omp


#endif /* CIRCLE_BOUND_H_ */
