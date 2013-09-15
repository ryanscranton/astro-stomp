/*
 * annulus_bound.h
 *
 *  Created on: Jul 5, 2012
 *      Author: scranton
 */

#ifndef ANNULUS_BOUND_H_
#define ANNULUS_BOUND_H_

#include "bound_interface.h"
#include "circle_bound.h"
#include "point.h"

namespace s2omp {

class annulus_bound: public bound_interface {
public:
  annulus_bound();
  annulus_bound(const point& axis, double inner_height, double outer_height);
  virtual ~annulus_bound();

  static annulus_bound* from_radii(const point& axis,
      double inner_theta_degrees, double outer_theta_degrees);
  static annulus_bound* from_heights(const point& axis, double inner_height,
      double outer_height);
  static annulus_bound* from_angular_bin(const point& axis,
      const angular_bin& bin);

  // Getters for the basic bound parameters.
  inline point axis() const {
    return axis_;
  }
  inline double inner_radius() const {
    return 2.0 * asin(sqrt(0.5 * inner_height_)) * RAD_TO_DEG;
  }
  inline double outer_radius() const {
    return 2.0 * asin(sqrt(0.5 * outer_height_)) * RAD_TO_DEG;
  }
  inline double inner_height() const {
    return inner_height_;
  }
  inline double outer_height() const {
    return outer_height_;
  }

  // Change the central axis
  inline void set_axis(const point& p) {
    axis_ = p;
    update_caps();
  }

  inline bool is_valid() const {
    return double_le(inner_height_, outer_height_) && double_le(0.0,
        inner_height_) && double_le(outer_height_, 2.0);
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

  virtual point get_center() const {
    return axis_;
  }
  virtual circle_bound get_bound() const {
    return circle_bound(axis_, outer_height_);
  }

private:
  inline void update_caps();
  inline bool inner_contains(const point& p) const;
  inline bool inner_contains(const pixel& pix) const;
  inline bool inner_may_intersect(const pixel& pix) const;
  inline bool outer_contains(const point& p) const;
  inline bool outer_contains(const pixel& pix) const;
  inline bool outer_may_intersect(const pixel& pix) const;

  point axis_;
  double inner_height_, outer_height_;
  S2Cap inner_cap_, outer_cap_;
};

void annulus_bound::update_caps() {
  inner_cap_ = S2Cap::FromAxisHeight(axis_.s2point(), inner_height_);
  outer_cap_ = S2Cap::FromAxisHeight(axis_.s2point(), outer_height_);
}

bool annulus_bound::inner_contains(const point& p) const {
  return double_le((axis_.s2point() - p.s2point()).Norm2(),
      2.0 * inner_height_);
}

bool annulus_bound::inner_contains(const pixel& pix) const {
  return inner_cap_.Contains(pix.get_cell());
}

bool annulus_bound::inner_may_intersect(const pixel& pix) const {
  return inner_cap_.MayIntersect(pix.get_cell());
}

bool annulus_bound::outer_contains(const point& p) const {
  return double_le((axis_.s2point() - p.s2point()).Norm2(),
      2.0 * outer_height_);
}

bool annulus_bound::outer_contains(const pixel& pix) const {
  return outer_cap_.Contains(pix.get_cell());
}

bool annulus_bound::outer_may_intersect(const pixel& pix) const {
  return outer_cap_.MayIntersect(pix.get_cell());
}

} // end namespace s2omp

#endif /* ANNULUS_BOUND_H_ */
