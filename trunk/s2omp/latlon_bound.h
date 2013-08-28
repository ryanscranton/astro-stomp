/*
 * latlon_bound.h
 *
 *  Created on: Jul 5, 2012
 *      Author: scranton
 */

#ifndef LATLON_BOUND_H_
#define LATLON_BOUND_H_

#include <stdint.h>
#include <math.h>
#include <vector>
#include <algorithm>

#include <s2/s2.h>
#include <s2/s2latlngrect.h>

#include "core.h"
#include "point.h"
#include "bound_interface.h"

namespace s2omp {

class latlon_bound: public bound_interface {
public:
  latlon_bound();
  latlon_bound(const point& hi, const point& lo);
  virtual ~latlon_bound();

  // Special case static constructors for bands of constant RA or DEC.
  static latlon_bound from_ra_bounds(double ra_min_degrees,
      double ra_max_degrees);
  static latlon_bound from_dec_bounds(double dec_min_degrees,
      double dec_max_degrees);

  // Expand the current latlon_bound to include the input point.
  void add_point(const point& p);

  // Accessor for the bound vertices.
  point vertex(int k) const;

  // Returns true if the latitudinal and longitudinal bounds are valid.
  bool is_valid() const;

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

private:
  latlon_bound(S2LatLngRect bound);
  S2LatLngRect latlon_;
};

} // end namespace s2omp

#endif /* LATLON_BOUND_H_ */
