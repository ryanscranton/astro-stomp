/*
 * polygon_bound.h
 *
 *  Created on: Jul 5, 2012
 *      Author: scranton
 */

#ifndef POLYGON_BOUND_H_
#define POLYGON_BOUND_H_

namespace s2omp {

class polygon_bound : public geometric_bound {
public:
  virtual ~polygon_bound();

  static polygon_bound* from_points(const point_vector& points);

  bool add_loop(const point_vector& points);

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
