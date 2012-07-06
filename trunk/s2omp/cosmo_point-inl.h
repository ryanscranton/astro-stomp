// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the classes used to store point-like data on the
// sphere.  This can be simply a location on the sky (the AngularCoordinate
// class) or a location along with additional information about the object at
// that position (WeightedAngularCoordinate).  CosmoCoordinate extends the
// latter functionality to also allow for treatment of cosmological coordinates.

#ifndef COSMO_POINT_H_
#define COSMO_POINT_H_

#include "point.h"
#include "util.h"

namespace s2omp {

class cosmo_point: public point {
public:
	cosmo_point(double x, double y, double z, double weight, double redshift);

	static cosmo_point* from_radec_deg(double ra_deg, double dec_deg,
			double redshift);
	static cosmo_point* from_radec_deg(double ra_deg, double dec_deg,
			double weight, double redshift);
	static cosmo_point* from_radec_rad(double ra_rad, double dec_rad,
			double redshift);
	static cosmo_point* from_radec_rad(double ra_rad, double dec_rad,
			double weight, double redshift);
	static cosmo_point* from_point(const point& p, double redshift);

	double redshift() const;
	void set_redshift(double redshift);

	// Since we have a redshift attached to the location, we can now attach some
	// functionality related to the 3-D coordinates.  Start with some conversions
	// between angular distance between an input point and our coordinate and
	// the projected radius at our coordinate's redshift.  Since we're getting
	// this data out of the cosmology class, we follow that class's convention
	// of giving distance in comoving Mpc/h
	double projected_radius(const point& p) const;
	double projected_radius(const point* p) const;

	// We can also offer variations on the dot product that takes
	// CosmoCoordinates and do their calculations with the full 3-D vectors
	double dot(const cosmo_point& p) const;
	double dot(const cosmo_point* p) const;

	// Some methods for accessing the distance values implied by our redshift.
	double comoving_distance() const;
	double angular_diameter_distance() const;
	double luminosity_distance() const;

private:
	double redshift_;
};

inline cosmo_point::cosmo_point(double x, double y, double z, double weight,
		double redshift) {
	set_xyz(x, y, z);
	set_weight(weight);
	redshift_ = redshift;
}

inline cosmo_point* cosmo_point::from_radec_deg(double ra_deg, double dec_deg,
		double redshift) {
	return from_radec_rad(ra_deg * DEG_TO_RAD, dec_deg * DEG_TO_RAD, 1.0,
			redshift);
}

inline cosmo_point* cosmo_point::from_radec_deg(double ra_deg, double dec_deg,
		double weight, double redshift) {
	return from_radec_rad(ra_deg * DEG_TO_RAD, dec_deg * DEG_TO_RAD, weight,
			redshift);
}

inline cosmo_point* cosmo_point::from_radec_rad(double ra_rad, double dec_rad,
		double redshift) {
	return from_radec_rad(ra_rad, dec_rad, 1.0, redshift);
}

inline cosmo_point* cosmo_point::from_radec_rad(double ra_rad, double dec_rad,
		double weight, double redshift) {
	cosmo_point* p = new cosmo_point(0.0, 0.0, 0.0, weight, redshift);
	p->set_latlon_rad(dec_rad, ra_rad, EQUATORIAL);
	return p;
}

inline cosmo_point* cosmo_point::from_point(const point& p, double redshift) {
	return new cosmo_point(p.unit_sphere_x(), p.unit_sphere_y(),
			p.unit_sphere_z(), p.weight(), redshift);
}

inline double cosmo_point::redshift() {
	return redshift_;
}

inline void cosmo_point::set_redshift(double redshift) {
	redshift_ = redshift;
}

inline double cosmo_point::projected_radius(const point& p) {
	return cosmology::projected_distance(redshift_, angular_distance(p));
}

inline double cosmo_point::projected_radius(const point* p) {
	return cosmology::projected_distance(redshift_, angular_distance(p));
}

inline double cosmo_point::dot(const cosmo_point& p) {
	return cosmology::comoving_distance(redshift_) * (unit_sphere_x()
			* p.unit_sphere_x() + unit_sphere_y() * p.unit_sphere_y()
			+ unit_sphere_z() * p.unit_sphere_z());
}

inline double cosmo_point::dot(const cosmo_point* p) {
	return cosmology::comoving_distance(redshift_) * (unit_sphere_x()
			* p->unit_sphere_x() + unit_sphere_y() * p->unit_sphere_y()
			+ unit_sphere_z() * p->unit_sphere_z());
}

inline double cosmo_point::comoving_distance() {
	return cosmology::comoving_distance(redshift_);
}

inline double cosmo_point::angular_diameter_distance() {
	return cosmology::angular_diameter_distance(redshift_);
}

inline double cosmo_point::luminosity_distance() {
	return cosmology::luminosity_distance(redshift_);
}

} // end namespace s2omp

#endif /* COSMO_POINT_H_ */
