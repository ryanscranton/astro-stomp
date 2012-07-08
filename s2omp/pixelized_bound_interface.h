// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the abstract BaseMap class that serves as the
// basis for all of the *Map objects.  BaseMap sets out the basic functionality
// that all of the *Map classes need to describe a given region on the sky and
// do some basic internal maintenance.  Additionally, BaseMap provides a
// common set of methods for dividing that area up into nearly equal-area,
// similarly-shaped regions.  This functionality is the basis for calculating
// jack-knife errors for our various statistical analyses.

#ifndef PIXELIZED_BOUND_INTERFACE_H_
#define PIXELIZED_BOUND_INTERFACE_H_

typedef std::map<const uint64, int16_t> region_dict;
typedef region_dict::iterator region_iterator;
typedef std::pair<region_iterator, region_iterator> region_pair;
typedef std::map<const int16_t, double> region_area;

namespace s2omp {
class region_map {
	// This class provides the functionality for dividing the area subtended by a
	// BaseMap-derived object into roughly equal-area, equal-sized regions.  The
	// class is not intended to be instantiated outside of the BaseMap class.

public:
	region_map();
	virtual ~region_map();

	uint16_t initialize(uint16_t n_region, pixelized_bound* bound);
	uint16_t
			initialize(uint16_t n_region, uint32_t level, pixelized_bound* bound);

	bool initialize(pixelized_bound& reference_bound, pixelized_bound* bound);

	int16_t find_region(const point& p);

	int16_t find_region(const pixel& pix);

	void clear();

	void region_covering(int16_t region, pixel_vector* pixels);

	// Given a region index, return the area associated with that region.
	double region_area(int16_t region);

	// Some getter methods to describe the state of the RegionMap.
	inline uint16_t n_region() const {
		return n_region_;
	}
	inline uint32_t level() const {
		return level_;
	}
	inline bool is_initialized() const {
		return !region_map_.empty();
	}

	// Return iterators for the set of RegionMap objects.
	inline region_iterator begin() {
		return region_map_.begin();
	}
	inline region_iterator end() {
		return region_map_.end();
	}

private:
	int find_regionation_level(analytic_bound* bound, uint16_t n_region);

	void regionate(analytic_bound* bound, uint16_t n_region, int level);

	bool verify_regionation(uint16_t n_region);

	region_dict region_map_;
	region_area region_area_;
	int level_;
	uint16_t n_region_;
};

class pixelized_bound_interface {
public:
	virtual ~pixelized_bound_interface();

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

	// These methods all act as wrappers for the RegionMapper object contained
	// in the class.  See that class for documentation.
	inline uint16_t initialize_regions(uint16_t n_regions) {
		return region_map_.initialize_regions(n_regions);
	}
	inline uint16_t initialize_regions(uint16_t n_regions, int region_level) {
		return region_map_.initialize_regions(n_regions, region_level);
	}
	inline bool initialize_regions(pixelized_bound* bound) {
		return region_map_.initialize_regions(bound);
	}
	inline int16_t find_region(const point& p) const {
		return region_map_.find_region(p);
	}
	inline int16_t find_region(const pixel& pix) const {
		return region_map_.find_region(pix);
	}
	inline void clear_regions() {
		return region_map_.clear();
	}
	inline void region_covering(int16_t region, pixel_vector* pixels) const {
		region_map_.region_covering(region, pixels);
	}
	inline double region_area(int16_t region) const {
		return region_map_.region_arae(region);
	}
	inline uint16_t n_region() const {
		return region_map_.n_region();
	}
	inline uint32_t region_level() const {
		return region_map_.level();
	}
	inline bool regions_initialized() const {
		return region_map_.is_initialized();
	}
	inline region_iterator region_begin() {
		return region_map_.begin();
	}
	inline region_terator region_end() {
		return region_map_.end();
	}

private:
	pixelized_bound_interface();
	region_map region_map_;
};

} // end namespace s2omp

#endif /* PIXELIZED_BOUND_INTERFACE_H_ */
