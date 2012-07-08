// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains a variant on the Pixel class.  The goal here is to
// encode some manner of scalar field (galaxy density on the sky, CMB
// temperature, etc.) where we need both the value and an associated noise on
// the value.  Likewise, as will be seen in the ScalarMap class, these pixels
// will generally be used to sample this scalar field uniformly over some
// region.  This is in contrast with the Map object where the goal is to
// accurately describe a region's geometry using pixels of various sizes.

#ifndef SCALAR_PIXEL_H_
#define SCALAR_PIXEL_H_

namespace s2omp {

class ScalarPixel;

typedef std::vector<scalar_pixel> scalar_vector;
typedef scalar_vector::iterator scalar_iterator;
typedef std::pair<scalar_iterator, scalar_iterator> scalar_pair;
typedef std::vector<scalar_pixel *> scalar_ptr_vector;
typedef scalar_ptr_vector::iterator scalar_ptr_iterator;

class scalar_pixel: public pixel {
public:
	explicit scalar_pixel(uint64 id);
	scalar_pixel(uint64 id, double intensity, double weight, uint32_t n_points);
	virtual ~scalar_pixel();

	static scalar_pixel from_pixel(pixel pix, double intensity, double weight,
			uint32_t n_points);
	static scalar_pixel from_point(point p, int level, double intensity,
			double weight, uint32_t n_points);

	inline double intensity() const {
		return intensity_;
	}
	inline double weight() const {
		return weight_;
	}
	inline uint32_t n_points() const {
		return n_points_;
	}
	inline double unit_sphere_x() const {
		return unit_sphere_x_;
	}
	inline double unit_sphere_y() const {
		return unit_sphere_y_;
	}
	inline double unit_sphere_z() const {
		return unit_sphere_z_;
	}

	inline void set_intensity(double intensity) :
		intensity_(intensity) {
	}
	inline void set_weight(double weight) :
		weight_(weight) {
	}
	inline void set_n_points(uint32_t n_points) :
		n_points_(n_points) {
	}

	double mean_intensity() const;
	void add_to_intensity(const double intensity, const uint32_t n_point);

	void convert_to_overdensity(double expected_intensity);

	void convert_to_fractional_overdensity(double expected_intensity);

	// And two complementary methods to take us from over-densities back to raw
	// intensities.
	void convert_from_overdensity(double expected_intensity);
	void convert_from_fractional_overdensity(double expected_intensity);

	// Finally, a method to tell us whether we're in over-density mode or not.
	bool is_overdensity() const {
		return is_overdensity_;
	}

private:
	scalar_pixel();
	// An internal method that we'll use when we calculate correlation
	// functions.
	void _within_annulus(angular_bin& theta, scalar_vector& pix);

	double intensity_, weight_;
	double unit_sphere_x_, unit_sphere_y_, unit_sphere_z_;
	uint32_t n_points_;
	bool is_overdensity_;
};

} // end namespace s2omp

#endif /* SCALAR_PIXEL_H_ */
