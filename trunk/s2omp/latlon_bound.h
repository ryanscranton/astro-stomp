/*
 * latlon_bound.h
 *
 *  Created on: Jul 5, 2012
 *      Author: scranton
 */

#ifndef LATLON_BOUND_H_
#define LATLON_BOUND_H_

namespace s2omp {

class latlon_bound : public geometric_bound {
public:
  virtual ~latlon_bound();

  static latlon_bound* from_points(const point& hi, const point& lo);

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
  latlon_bound();
  S2LatLonRect latlon_;
};

} // end namespace s2omp

#endif /* LATLON_BOUND_H_ */
