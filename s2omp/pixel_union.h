// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the Map class.  Maps are intended to describe an
// arbitrary area on the sky as precisely as possible, given the limits of
// pixel resolution, physical memory, etc.  Internally, this information is
// encoded using pixels at a range of resolutions, depending on how the
// boundaries of the area in question and the pixelization scheme interact.
// However, the goal of the class is to abstract away those details, allowing
// the user to treat Maps as a pure representative of spherical geometry.

#ifndef PIXEL_UNION_H_
#define PIXEL_UNION_H_

namespace s2omp {

class pixel_union: public pixelized_bound_interface {
public:
	pixel_union();
	virtual ~pixel_union();

	static pixel_union* from_covering(const pixel_vector& pixels);

	void init(const pixel_vector& pixels);

	int min_level() const;
	int max_level() const;

	void soften(int max_level);

	void combine(const pixel_union& u);
	void intersect(const pixel_union& u);
	void exclude(const pixel_union& u);

	void init_from_combination(const pixel_union& a, const pixel_union& b);
	void init_from_intersection(const pixel_union& a, const pixel_union& b);
	void init_from_exclusion(const pixel_union& a, const pixel_union& b);

	point get_random_point() const;
	void get_random_points(long n_points, pixel_vector* points) const;

	pixel_iterator begin();
	pixel_iterator end();

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
	pixel_vector pixels_;
	int min_level_, max_level_;
	double area_;
	bool initialized_;
};

} // end namespace s2omp

#endif /* PIXEL_UNION_H_ */
