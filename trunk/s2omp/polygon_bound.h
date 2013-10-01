/*
 * polygon_bound.h
 *
 *  Created on: Jul 5, 2012
 *      Author: scranton
 */

#ifndef POLYGON_BOUND_H_
#define POLYGON_BOUND_H_

#include <stdint.h>
#include <math.h>
#include <vector>
#include <algorithm>

#include <s2/s2.h>
#include <s2/s2polygon.h>

#include "core.h"
#include "point.h"
#include "bound_interface.h"

namespace s2omp {

class polygon_bound : public geometric_bound {
public:
  polygon_bound(const point_vector& points);
  virtual ~polygon_bound();

  static polygon_bound* from_points(const point_vector& points);

  bool add_loop(const point_vector& points);

  inline bool is_valid() const {
    return polygon_.IsValid();
  }

  bool intersects(const polygon_bound& bound) const;

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
  polygon_bound();
  S2Polygon polygon_;
};

} // end namespace s2omp

#endif /* POLYGON_BOUND_H_ */
