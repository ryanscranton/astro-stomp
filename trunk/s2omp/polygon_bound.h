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

class polygon_bound : public bound_interface {
public:
  polygon_bound(const point_vector& points);
  virtual ~polygon_bound();

  static polygon_bound* from_points(const point_vector& points);

  bool add_loop(const point_vector& points);

  // We require a non-empty polygon for our validity check, unlike S2Polygon.
  inline bool is_valid() const {
    return polygon_->num_vertices() > 0 && polygon_->IsValid();
  }

  inline long num_vertices() const {
    return polygon_->num_vertices();
  }

  inline long num_loops() const {
    return polygon_->num_loops();
  }

  bool intersects(const polygon_bound& bound) const;

  static S2Loop* point_vector_to_s2loop(const point_vector& points);

  // API from bound_interface.h
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

protected:
  inline S2Polygon* get_s2polygon() const;

private:
  polygon_bound();
  S2Polygon* polygon_;
};

S2Polygon* polygon_bound::get_s2polygon() const {
  return polygon_;
}

} // end namespace s2omp

#endif /* POLYGON_BOUND_H_ */
