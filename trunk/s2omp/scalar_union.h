// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the ScalarMap class.  Unlike Maps, the primary
// goal here is to encode a scalar field over some area of the sky.  As such,
// we sacrifice some degree of precision in describing the exact area of the
// field and we use a uniform sampling of the field across the area in
// question.  This makes the class ideal for calculating angular correlation
// functions on the encoded field.


#ifndef SCALAR_UNION_H_
#define SCALAR_UNION_H_

namespace s2omp {

class scalar_union: public pixelized_bound_interface {
public:
	enum ScalarType {
		SCALAR_FIELD, DENSITY_FIELD, SAMPLED_FIELD
	};

	explicit scalar_union(ScalarType t);
	~scalar_union();

	static scalar_union* from_bound(const geometric_bound& bound, int level,
			ScalarType t);
	static scalar_union* from_pixel_union(const pixel_union& p, int level,
			ScalarType t);
	static scalar_union* from_scalar_pixels(const scalar_vector& pixels,
			ScalarType t);
	static scalar_union* from_scalar_union(const scalar_union& s, int level);

	bool init(const geometric_bound& bound, int level, ScalarType t);
	bool init(const pixel_union& p, int level, ScalarType t);
	bool init(const scalar_vector& pixels, ScalarType t);

	bool add_point(point& p);

	scalar_pixel resample(const pixel& p) const;

	double find_intensity(const pixel& p) const;
	double find_density(const pixel& p) const;
	double find_point_density(const pixel& p) const;

	double find_local_area(annulus_bound& bound);
	double find_local_intensity(annulus_bound& bound);
	double find_local_density(annulus_bound& bound);
	double find_local_point_density(annulus_bound& bound);

	void calculate_mean_intensity();
	void convert_to_over_density();
	void convert_from_over_density();

	void auto_correlate(theta_iterator theta_iter);
	void auto_correlate(angular_orrelation* wtheta);

	void auto_correlate_with_regions(angular_correlation* wtheta);
	void auto_correlate_with_regions(theta_iterator theta_iter);

	void cross_correlate(const scalar_union* s, angular_correlation* wtheta);
	void cross_correlate(const scalar_union* s, theta_iterator theta_iter);
	void cross_correlate_with_regions(scalar_union* s,
			angular_correlation* wtheta);
	void cross_correlate_with_regions(scalar_union& s, theta_iterator theta_iter);

	int level() const;
	ScalarType type() const;
  double intensity() const;
  int n_points() const;
  double density() const;
  double point_density() const;
  double mean_intensity() const;
  bool is_over_density();
  scalar_iterator begin();
  scalar_iterator end();

  // API from pixelized_bound_interface.h
	virtual bool is_empty() const;
	virtual long size() const;
	virtual void clear() const;
	virtual void area() const;

	virtual bool contains(const point& p) const;
	virtual bool contains(const pixel& pix) const;

	virtual double contained_area(const pixel& pix) const;
	virtual bool may_intersect(const pixel& pix) const;

	virtual void covering(pixel_vector* pixels) const;
	virtual void covering(int max_pixels, pixel_vector* pixels) const;
	virtual void simple_covering(int level, pixel_vector* pixels) const;

	virtual circle_bound* get_bound() const;

private:
	scalar_union();

	scalar_vector pixels_;
  double area_, mean_intensity_, unmasked_fraction_minimum_, total_intensity_;
	ScalarType type_;
	int level_;
	uint32_t total_points_;
	bool initialized_, converted_to_overdensity_, calculated_mean_intensity_;
};

} // end namespace s2omp

#endif /* SCALAR_UNION_H_ */
