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

class pixel;  // class declaration in pixel.h
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
    Equatorial,
    Galactic,
    Geocentric,
    Heliocentric
  };

 public:
  static point* from_latlon_degrees(double lat_deg, double lon_deg, Sphere s);
  static point* from_latlon_radians(double lat_rad, double lon_rad, Sphere s);

  static point* from_radec_degrees(double ra_deg, double dec_deg);
  static point* from_radec_radians(double ra_rad, double dec_deg);

  double lat_deg(Sphere s);
  double lon_deg(Sphere s);
  double lat_rad(Sphere s);
  double lon_rad(Sphere s);

  double ra_deg();
  double dec_deg();
  double ra_rad();
  double dec_rad();

  double unit_sphere_x();
  double unit_sphere_y();
  double unit_sphere_z();

  double angular_distance(point& p);
  double angular_distance(point* p);

  double dot(point& p);
  double dot(point* p);
  static double dot(point& a, point& b);

  point* cross(point& p);
  point* cross(point* p);
  static point* cross(point& a, point& b);

  point* great_circle(point& p, Sphere s);
  static point* great_circle(point& a, point& b, Sphere s);

  double position_angle(point& p, Sphere s);

  void rotate_about(point& axis, double rotation_angle_degrees, Sphere s);
  static point* rotate_about(point& p, point& axis, double rotation_angle, Sphere s);

  pixel* to_pixel();
  pixel* to_pixel(int level);

 private:
  point();
  virtual ~point();

  point(double x, double y, double z);
  double cos_position_angle(point& p, Sphere s);
  double sin_position_angle(point& p, Sphere s);
  void rotate_about(point& axis, double rotation_angle_degrees,
      Sphere s, double& unit_sphere_x, double& unit_sphere_y,
                double& unit_sphere_z);

  S2::SPoint point_;

};

} // end namespace s2omp
#endif /* POINT_H_ */
