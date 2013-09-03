#include <gtest/gtest.h>

#include "pixel.h"
#include "point.h"

TEST(point, TestPointDefaultConstructor) {
  // The default constructor without any arguments should produce an invalid
  // point.
  s2omp::point p;
  ASSERT_FALSE(p.is_normalized());

  // Create a point with the x-y-z constructor
  p = s2omp::point(0.0, 0.0, 1.0, 1.0);
  ASSERT_TRUE(p.is_normalized());
  ASSERT_DOUBLE_EQ(p.x(), 0.0);
  ASSERT_DOUBLE_EQ(p.y(), 0.0);
  ASSERT_DOUBLE_EQ(p.z(), 1.0);
  ASSERT_DOUBLE_EQ(p.weight(), 1.0);

  // Verify that the constructor properly normalizes a point even if the input
  // x-y-z is not unit normal.
  p = s2omp::point(0.0, 0.0, 2.0, 1.0);
  ASSERT_TRUE(p.is_normalized());
  ASSERT_DOUBLE_EQ(p.x(), 0.0);
  ASSERT_DOUBLE_EQ(p.y(), 0.0);
  ASSERT_DOUBLE_EQ(p.z(), 1.0);

  p = s2omp::point(1.0, 1.0, 1.0, 1.0);
  ASSERT_TRUE(p.is_normalized());
}

TEST(point, TestPointCoordinateConstructors) {
  // Test the various constructors that take input coordinates in the various
  // coordinate systems.
  double ra = 0.0;
  double dec = 90.0;
  double expected_precision = 1.0e-6;

  s2omp::point p = s2omp::point::from_radec_deg(ra, dec, 1.0);
  ASSERT_TRUE(p.is_normalized());
  ASSERT_DOUBLE_EQ(p.ra_deg(), ra);
  ASSERT_DOUBLE_EQ(p.dec_deg(), dec);
  ASSERT_DOUBLE_EQ(p.lat_deg(s2omp::point::EQUATORIAL), dec);
  ASSERT_DOUBLE_EQ(p.lon_deg(s2omp::point::EQUATORIAL), ra);
  ASSERT_DOUBLE_EQ(p.weight(), 1.0);
  ASSERT_NEAR(p.x(), 0.0, expected_precision);
  ASSERT_NEAR(p.y(), 0.0, expected_precision);
  ASSERT_NEAR(p.z(), 1.0, expected_precision);
  ASSERT_NEAR(p.unit_sphere_x(), 0.0, expected_precision);
  ASSERT_NEAR(p.unit_sphere_y(), 0.0, expected_precision);
  ASSERT_NEAR(p.unit_sphere_z(), 1.0, expected_precision);

  double glat = 0.0;
  double glon = 180.0;
  p = s2omp::point::from_latlon_deg(glat, glon,
      s2omp::point::Coordinate::GALACTIC);
  ASSERT_TRUE(p.is_normalized());
  ASSERT_NEAR(p.lat_deg(s2omp::point::GALACTIC), glat, expected_precision);
  ASSERT_NEAR(p.lon_deg(s2omp::point::GALACTIC), glon, expected_precision);
  ASSERT_NEAR(p.unit_sphere_x(s2omp::point::GALACTIC), -1.0,expected_precision);
  ASSERT_NEAR(p.unit_sphere_y(s2omp::point::GALACTIC), 0.0, expected_precision);
  ASSERT_NEAR(p.unit_sphere_z(s2omp::point::GALACTIC), 0.0, expected_precision);

  double elat = 0.0;
  double elon = 0.0;
  p = s2omp::point::from_latlon_deg(elat, elon, s2omp::point::ECLIPTIC);
  ASSERT_TRUE(p.is_normalized());
  ASSERT_NEAR(p.lat_deg(s2omp::point::ECLIPTIC), elat, expected_precision);
  ASSERT_NEAR(p.lon_deg(s2omp::point::ECLIPTIC), elon, expected_precision);
  ASSERT_NEAR(p.unit_sphere_x(s2omp::point::ECLIPTIC), 1.0, expected_precision);
  ASSERT_NEAR(p.unit_sphere_y(s2omp::point::ECLIPTIC), 0.0, expected_precision);
  ASSERT_NEAR(p.unit_sphere_z(s2omp::point::ECLIPTIC), 0.0, expected_precision);
}

TEST(point, TestPointCoordinateConversion) {
  // Test conversion between the various coordinate systems.

  double ra_rad = 45.0 * s2omp::DEG_TO_RAD;
  double dec_rad = 45.0 * s2omp::DEG_TO_RAD;

  // Expected coordinates from http://ned.ipac.caltech.edu/forms/calculator.html
  double expected_glon_rad = 145.56300265 * s2omp::DEG_TO_RAD;
  double expected_glat_rad = -12.14813638 * s2omp::DEG_TO_RAD;
  double expected_elon_rad = 55.95449205 * s2omp::DEG_TO_RAD;
  double expected_elat_rad = 26.73529289 * s2omp::DEG_TO_RAD;

  double glon_rad, glat_rad;
  double expected_precision = 1.0e-5;
  s2omp::point::equatorial_to_galactic(ra_rad, dec_rad, &glon_rad, &glat_rad);
  ASSERT_NEAR(expected_glon_rad, glon_rad, expected_precision);
  ASSERT_NEAR(expected_glat_rad, glat_rad, expected_precision);

  double elon_rad, elat_rad;
  s2omp::point::equatorial_to_ecliptic(ra_rad, dec_rad, &elon_rad, &elat_rad);
  ASSERT_NEAR(expected_elon_rad, elon_rad, expected_precision);
  ASSERT_NEAR(expected_elat_rad, elat_rad, expected_precision);

  double converted_ra_rad, converted_dec_rad;
  s2omp::point::galactic_to_equatorial(expected_glon_rad, expected_glat_rad,
      &converted_ra_rad, &converted_dec_rad);
  ASSERT_NEAR(ra_rad, converted_ra_rad, expected_precision);
  ASSERT_NEAR(dec_rad, converted_dec_rad, expected_precision);

  s2omp::point::ecliptic_to_equatorial(expected_elon_rad, expected_elat_rad,
      &converted_ra_rad, &converted_dec_rad);
  ASSERT_NEAR(ra_rad, converted_ra_rad, expected_precision);
  ASSERT_NEAR(dec_rad, converted_dec_rad, expected_precision);
}

TEST(point, TestPointDotProducts) {
  // Test routines for calculating dot products.

  s2omp::point p = s2omp::point::from_radec_deg(0.0, 0.0);
  s2omp::point orth = s2omp::point::from_radec_deg(0.0, 90.0);
  double expected_precision = 1.0e-10;

  ASSERT_NEAR(p.dot(p), 1.0, expected_precision);
  ASSERT_NEAR(s2omp::point::dot(p, p), 1.0, expected_precision);
  ASSERT_NEAR(p.dot(orth), 0.0, expected_precision);
}

TEST(point, TestPointCrossProducts) {
  // Test routines for calculating cross products

  s2omp::point p = s2omp::point::from_radec_deg(0.0, 0.0);
  s2omp::point orth = s2omp::point::from_radec_deg(0.0, 90.0);
  double expected_precision = 1.0e-10;

  s2omp::point p_cross_orth = p.cross(orth);
  s2omp::point p_cross_orth_static = s2omp::point::cross(p, orth);
  s2omp::point orth_cross_p_static = s2omp::point::cross(orth, p);
  ASSERT_NEAR(p.dot(p_cross_orth), 0.0, expected_precision);
  ASSERT_NEAR(s2omp::point::dot(p, p_cross_orth_static), 0.0,
      expected_precision);
  ASSERT_NEAR(s2omp::point::dot(p, orth_cross_p_static), 0.0,
      expected_precision);
  ASSERT_NEAR(s2omp::point::dot(orth_cross_p_static, p_cross_orth_static), -1.0,
      expected_precision);
}

TEST(point, TestPointRotation) {
  // Test routines for rotating points along great circles.

  s2omp::point p = s2omp::point::from_radec_deg(0.0, 0.0);
  s2omp::point orth = s2omp::point::from_radec_deg(0.0, 90.0);
  double expected_precision = 1.0e-10;
  double theta_deg = 90.0;

  s2omp::point rotated_p = s2omp::point::rotate_about(p, orth, theta_deg,
      s2omp::point::EQUATORIAL);
  ASSERT_NEAR(p.dot(rotated_p), cos(theta_deg * s2omp::DEG_TO_RAD),
      expected_precision);
  ASSERT_NEAR(rotated_p.dot(orth), 0.0, expected_precision);

  theta_deg = 180.0;
  rotated_p = s2omp::point::rotate_about(p, orth, theta_deg,
      s2omp::point::EQUATORIAL);
  ASSERT_NEAR(p.dot(rotated_p), cos(theta_deg * s2omp::DEG_TO_RAD),
      expected_precision);
}

TEST(point, TestPointPositionAngles) {
  // Test routines for calculating position angles.

  s2omp::point p = s2omp::point::from_radec_deg(0.0, 0.0);
  s2omp::point north = s2omp::point::from_radec_deg(0.0, 10.0);
  s2omp::point east = s2omp::point::from_radec_deg(10.0, 0.0);
  double expected_precision = 1.0e-10;

  ASSERT_NEAR(p.position_angle(east, s2omp::point::EQUATORIAL),
      90.0, expected_precision);
  ASSERT_NEAR(s2omp::point::position_angle(p, east,
      s2omp::point::EQUATORIAL), 90.0, expected_precision);
  ASSERT_NEAR(p.position_angle(north, s2omp::point::EQUATORIAL),
      0.0, expected_precision);
}
