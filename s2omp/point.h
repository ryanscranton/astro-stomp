// Copyright 2012  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the classes used to store point-like data on the
// sphere.  This can be simply a location on the sky (the point
// class) or a location along with additional information about the object at
// that position (Weightedpoint).  CosmoCoordinate extends the
// latter functionality to also allow for treatment of cosmological coordinates.


#ifndef POINT_H_
#define POINT_H_

#include <math.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>

namespace s2omp {

class pixel; // class declaration in pixel.h
class point;

typedef std::vector<point> point_vector;
typedef point_vector::iterator point_iterator;
typedef std::vector<point *> point_ptr_vector;
typedef point_ptr_vector::iterator point_ptr_iterator;

class point {
	// Our generic class for handling angular positions.  The idea is that
	// locations on the celestial sphere should be abstract objects from which
	// you can draw whatever angular coordinate pair is necessary for a given
	// use case.  point's can be instantiated with a particular
	// coordinate system in mind or that can be set later on.
public:
	enum Sphere {
		EQUATORIAL, GALACTIC, GEOCENTRIC, HELIOCENTRIC
	};

public:
	point(double x, double y, double z, double weight);
	virtual ~point();

	static point from_latlon_degrees(double lat_deg, double lon_deg, Sphere s);
	static point from_latlon_degrees(double lat_deg, double lon_deg, Sphere s,
			double weight);
	static point from_latlon_radians(double lat_rad, double lon_rad, Sphere s);
	static point from_latlon_radians(double lat_rad, double lon_rad, Sphere s,
			double weight);

	static point from_radec_degrees(double ra_deg, double dec_deg);
	static point from_radec_degrees(double ra_deg, double dec_deg, double weight);
	static point from_radec_radians(double ra_rad, double dec_rad);
	static point from_radec_radians(double ra_rad, double dec_rad, double weight);

	void set_latlon_degrees(double lat_deg, double lon_deg, Sphere s);
	void set_latlon_radians(double lat_rad, double lon_rad, Sphere s);
	void set_xyz(double x, double y, double z);

	double lat_deg(Sphere s) const;
	double lon_deg(Sphere s) const;
	double lat_rad(Sphere s) const;
	double lon_rad(Sphere s) const;

	double ra_deg() const;
	double dec_deg() const;
	double ra_rad() const;
	double dec_rad() const;

	double unit_sphere_x() const;
	double unit_sphere_y() const;
	double unit_sphere_z() const;

	double weight() const;
	double set_weight(double weight);

	double angular_distance(const point& p) const;
	double angular_distance(const point* p) const;
	static double angular_distance(const point& a, const point& b);

	double dot(const point& p);
	double dot(const point* p);
	static double dot(const point& a, const point& b);

	point cross(const point& p);
	point cross(const point* p);
	static point cross(const point& a, const point& b);

	point great_circle(const point& p, Sphere s);
	static point great_circle(const point& a, const point& b, Sphere s);

	double position_angle(const point& p, Sphere s);
	static double position_angle(const point& center, const point& target,
			Sphere s);

	void rotate_about(point& axis, double rotation_angle_degrees, Sphere s);
	static point rotate_about(point& p, point& axis, double rotation_angle,
			Sphere s);

	pixel to_pixel();
	pixel to_pixel(int level);

private:
	point();

	double cos_position_angle(point& p, Sphere s);
	double sin_position_angle(point& p, Sphere s);
	void rotate_about(point& axis, double rotation_angle_degrees, Sphere s,
			double& unit_sphere_x, double& unit_sphere_y, double& unit_sphere_z);

	S2::SPoint point_;
	double weight_;
};

} // end namespace s2omp
#endif /* POINT_H_ */
