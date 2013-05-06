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
#include "point.h"
#include "bound_interface.h"
#include "MersenneTwister.h"

namespace s2omp {

class angular_bin;
class bound_interface;
class circle_bound;
class point;

typedef std::vector<circle_bound *> circle_ptr_vector;
typedef circle_ptr_vector::iterator circle_ptr_iterator;

class circle_bound : public bound_interface {
public:
  circle_bound();
  virtual ~circle_bound();

  // static circle_bound* from_angular_bin(const point& axis,
  //		const angular_bin& bin);
  static circle_bound* from_radius(const point& axis, double radius_degrees);
  static circle_bound* from_height(const point& axis, double height);

  inline point axis() {
    return axis_;
  }
  // inline void set_axis(const point& p);
  inline double radius() {
    return acos(1.0 - height_) * DEG_TO_RAD;
  }
  inline double height() {
    return height_;
  }

  // API from geometric_bound
  virtual bool is_empty() const;
  virtual long size() const;
  virtual void clear();
  virtual double area() const;

  virtual bool contains(const point& p) const;
  virtual bool contains(const pixel& pix) const;
  virtual double contained_area(const pixel& pix) const;

  virtual bool may_intersect(const pixel& pix) const;

  virtual point get_center() const;
  virtual circle_bound get_bound() const;

  virtual point get_random_point();
  virtual void get_random_points(long n_points, point_vector* points);
  virtual point get_weighted_random_point(const point_vector& points);
  virtual void get_weighted_random_points(
      long n_points, point_vector* points,
      const point_vector& input_points);

private:
  circle_bound(const point& axis, double height);

  bool intersects(const pixel& pix, const point_vector& vertices) const;
  circle_bound* get_complement() const;

  point axis_;
  double height_;

  // The variables below are for generating random points and will not be
  // initialized until initialize_random() is called.
  void initialize_random();
  point great_circle_norm_ ;
  double rotate_;
  bool initialized_random_;
  MTRand mtrand;
};

} // end namespace s2omp


#endif /* CIRCLE_BOUND_H_ */
