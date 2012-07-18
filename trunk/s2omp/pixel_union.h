// Copyright 2012  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the pixel_union class.  A pixel_union is analogous
// to the S2::CellUnion or Stomp::Map classes, where we describe a region of
// the sky with a collection of pixels.  Most of the functionality of the
// Map class is replicated here, with the exception of the fact that
// pixel_union objects do not encode scalar weights.

#ifndef PIXEL_UNION_H_
#define PIXEL_UNION_H_

#include <stdint.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include "core.h"
#include "point.h"
#include "pixel.h"
#include "pixelized_bound_interface.h"

namespace s2omp {

class pixel;
class point;

typedef std::map<const uint64, unit16_t> region_dict;
typedef pixel_map::iterator region_dict_iterator;

class pixel_union: public pixelized_bound_interface {
public:
	// Like the S2CellUnion and Stomp::Map classes, the pixels that comprise
	// a pixel_union need to be properly normalized to make the remainder of
	// the class work.  Since that can be time-consuming, we defer that
	// operation to the init() function.
	pixel_union();
	virtual ~pixel_union();

	// Alternatively, one can use this static method to return a pointer
	// to a normalized pixel_union from a vector of pixels.
	inline static pixel_union* from_covering(const pixel_vector& pixels) {
		pixel_union* p_union = new pixel_union();
		p_union->init(pixels);
		return p_union;
	}

	// Use the input vector of pixels to initialize this pixel_union.  These
	// pixels will be sorted by index and child pixels will be combined into
	// parent pixels where possible.
	void init(const pixel_vector& pixels);
	void init_swap(pixel_vector* pixels);

	// In some cases, we may wish to soften the edges of our pixel_union (for
	// instance if the current union is formed from combining multiple high
	// resolution unions producing a union that is unweidly).  In this case,
	// we combine all pixels smaller than the input level, producing a parent
	// pixel at max_level if more than half its area was in the original union.
	void soften(int max_level);

	// As with Stomp::Map, we have multiple methods for doing logical operations
	// on pixel_unions: union, intersection and exclusion.  In all cases, the
	// contents of the current pixel_union are replaced by the results of the
	// operation, which may lead to an empty pixel_union.
	void combine(const pixel_union& u);
	void intersect(const pixel_union& u);
	void exclude(const pixel_union& u);

	// Alternatively, we can initialize this pixel_union from the results of
	// doing a union, intersection or exclusion between two other pixel_unions.
	void init_from_combination(const pixel_union& a, const pixel_union& b);
	void init_from_intersection(const pixel_union& a, const pixel_union& b);
	void init_from_exclusion(const pixel_union& a, const pixel_union& b);

	// Return a Poisson-random point (or multiple such points) from the area
	// covered by this pixel_union
	point get_random_point() const;
	void get_random_points(long n_points, pixel_vector* points) const;

	inline int min_level() const {return min_level_;}
	inline int max_level() const {return max_level_;}

	inline pixel_iterator begin() {return pixels_.begin();}
	inline pixel_iterator end() {return pixels_.end();}

	// API from pixelized_bound_interface.h
	inline virtual bool is_empty() const {return !pixels_.empty();}
	inline virtual long size() const {return pixels_.size();}
	virtual void clear() const;
	inline virtual double area() const {return area_;}
	virtual double aveage_area() const;

	virtual bool contains(const point& p) const;
	virtual bool contains(const pixel& pix) const;
	virtual bool intersects(const pixel& pix) const;

	virtual double contained_area(const pixel& pix) const;

	virtual void covering(pixel_vector* pixels) const;
	virtual void covering(int max_pixels, pixel_vector* pixels) const;
	virtual void simple_covering(int level, pixel_vector* pixels) const;

	inline virtual circle_bound get_bound() const {
	  if (!initialized_bound_) initialize_bound();
	  return bound_;
	}

private:
	void generate_basic_covering(int level);
	void initialize_bound();

	void pixel_intersection(const pixel& pix, const pixel_union& u,
	    pixel_vector* pixels);
	void pixel_exclusion(const pixel& pix, const pixel_union& u,
	    pixel_vector* pixels);
	// For automatic regionation is is useful to know what the average area of
	// a pixel is in our union.

	pixel_vector pixels_;
	pixel_vector covering_;
	circle_bound bound_;
	int min_level_, max_level_;
	double area_;
	bool initialized_, initialized_bound_;
	region_map region_map_;
};

inline static pixel_union* pixel_union::from_covering(
		const pixel_vector& pixels) {
	pixel_union* p = new pixel_union();
	p->init(pixels);
	return p;
}

} // end namespace s2omp

#endif /* PIXEL_UNION_H_ */
