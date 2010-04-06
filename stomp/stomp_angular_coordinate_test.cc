#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>
#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "stomp_pixel.h"

void AngularCoordinateBasicTests() {
  // Ok, some basic routines.  We're going to declare an angular
  // position, convert it from Survey to Equatorial to Galactic coordinates
  // and make sure that it gives us the same pixel position regardless.
  std::cout << "\n";
  std::cout << "*************************************\n";
  std::cout << "*** AngularCoordinate Basic tests ***\n";
  std::cout << "*************************************\n";
  double lambda, eta, ra, dec, gal_lat, gal_lon;

  lambda = 10.0;
  eta = 10.0;

  Stomp::AngularCoordinate ang(lambda, eta, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);
  std::cout << "\tLambda: " << ang.Lambda() << ", Eta: " << ang.Eta() <<
    "; " << tmp_pix.HPixnum() << " " << tmp_pix.Superpixnum() << "\n";

  Stomp::AngularCoordinate::SurveyToEquatorial(lambda,eta,ra,dec);
  tmp_pix.SetPixnumFromAng(ang);
  std::cout << "\tRA: " << ang.RA() << " (" << ra << "), " <<
    "DEC: " << ang.DEC() <<  " (" << dec << "); " << tmp_pix.HPixnum() <<
    " " << tmp_pix.Superpixnum() << "\n";

  ang.SetSurveyCoordinates(lambda, eta);
  tmp_pix.SetPixnumFromAng(ang);
  Stomp::AngularCoordinate::SurveyToGalactic(lambda,eta,gal_lon,gal_lat);
  std::cout << "\tGalLat: " << ang.GalLat() << " (" << gal_lat << "), "
    "GalLon: " << ang.GalLon() <<  " (" << gal_lon << "); " <<
    tmp_pix.HPixnum() << " " << tmp_pix.Superpixnum() << "\n";

  Stomp::AngularCoordinate ang_a(0.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::AngularCoordinate ang_b(45.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::AngularCoordinate ang_c(60.0, 0.01, Stomp::AngularCoordinate::Survey);
  Stomp::AngularCoordinate ang_d(0.0, 60.0, Stomp::AngularCoordinate::Survey);
  std::cout << "\tPoint A: " << ang_a.Lambda() << ", " << ang_a.Eta() << "\n";
  std::cout << "\tPoint B: " << ang_b.Lambda() << ", " << ang_b.Eta() << "\n";
  std::cout << "\tPoint C: " << ang_c.Lambda() << ", " << ang_c.Eta() << "\n";
  std::cout << "\tPoint D: " << ang_d.Lambda() << ", " << ang_d.Eta() << "\n";
  std::cout << "\tA * B: " << ang_a.DotProduct(ang_b) << "\t\t\t B * A: " <<
    ang_b.DotProduct(ang_a) << "\n";
  Stomp::AngularCoordinate ang_n;
  ang_a.GreatCircle(ang_b, ang_n);
  Stomp::AngularCoordinate ang_m;
  ang_b.GreatCircle(ang_a, ang_m);

  std::cout << "\tA x B: " << ang_n.Lambda() << ", " << ang_n.Eta() <<
    "\t\t B x A: " << ang_m.Lambda() << ", " << ang_m.Eta() << "\n";
  std::cout << "\tDistance from AB to C: " <<
    asin(fabs(ang_n.DotProduct(ang_c)))/Stomp::DegToRad <<
    "; Distance from AB to D: " <<
    asin(fabs(ang_n.DotProduct(ang_d)))/Stomp::DegToRad << "\n";
  std::cout << "\tDistance from (A x B) to C: " <<
    ang_n.AngularDistance(ang_c) << "; Distance from (B x A) to D: " <<
    ang_m.AngularDistance(ang_d) << "\n";
}

void AngularCoordinatePositionAngleTests() {
  // Check to make sure that the position angle calculations are being done
  // correctly.
  std::cout << "\n";
  std::cout << "**********************************************\n";
  std::cout << "*** AngularCoordinate position angle check ***\n";
  std::cout << "**********************************************\n\n";
  double dec = 0.0;
  double ra = 0.0;
  double theta = 1.0;
  Stomp::AngularCoordinate ang(ra, dec, Stomp::AngularCoordinate::Equatorial);
  Stomp::AngularCoordinate ang2(ra+theta, dec,
				Stomp::AngularCoordinate::Equatorial);
  Stomp::Pixel pix(ang2, Stomp::MaxPixelResolution);
  std::cout << "(" << ang.RA() << ", " << ang.DEC() << ") -> angle (" <<
    ang2.RA() << ", " << ang2.DEC() << "): " << ang.PositionAngle(ang2) <<
    " degrees\n";
  std::cout << "(" << ang2.RA() << ", " << ang2.DEC() << ") -> angle (" <<
    ang.RA() << ", " << ang.DEC() << "): " << ang2.PositionAngle(ang) <<
    " degrees\n";
  std::cout << "(" << ang.RA() << ", " << ang.DEC() << ") -> pixel (" <<
    pix.RA() << ", " << pix.DEC() << "): " << ang.PositionAngle(pix) <<
    " degrees\n\n";

  ang2.SetEquatorialCoordinates(ra, dec+theta);
  pix.SetPixnumFromAng(ang2);
  std::cout << "(" << ang.RA() << ", " << ang.DEC() << ") -> angle (" <<
    ang2.RA() << ", " << ang2.DEC() << "): " << ang.PositionAngle(ang2) <<
    " degrees\n";
  std::cout << "(" << ang2.RA() << ", " << ang2.DEC() << ") -> angle (" <<
    ang.RA() << ", " << ang.DEC() << "): " << ang2.PositionAngle(ang) <<
    " degrees\n";
  std::cout << "(" << ang.RA() << ", " << ang.DEC() << ") -> pixel (" <<
    pix.RA() << ", " << pix.DEC() << "): " << ang.PositionAngle(pix) <<
    " degrees\n\n";

  ang2.SetEquatorialCoordinates(ra+theta, dec+theta);
  pix.SetPixnumFromAng(ang2);
  std::cout << "(" << ang.RA() << ", " << ang.DEC() << ") -> angle (" <<
    ang2.RA() << ", " << ang2.DEC() << "): " << ang.PositionAngle(ang2) <<
    " degrees\n";
  std::cout << "(" << ang2.RA() << ", " << ang2.DEC() << ") -> angle (" <<
    ang.RA() << ", " << ang.DEC() << "): " << ang2.PositionAngle(ang) <<
    " degrees\n";
  std::cout << "(" << ang.RA() << ", " << ang.DEC() << ") -> pixel (" <<
    pix.RA() << ", " << pix.DEC() << "): " << ang.PositionAngle(pix) <<
    " degrees\n\n";
}

void AngularCoordinateRotationTests() {
  // Check to make sure that the rotation calculations are being done correctly.
  std::cout << "\n";
  std::cout << "****************************************\n";
  std::cout << "*** AngularCoordinate rotation check ***\n";
  std::cout << "****************************************\n\n";
  double theta = 0.0;
  double phi = 0.0;
  double radius = 10.0;
  std::cout << "Starting in Equatorial coordinates.\n" <<
    "\tRotating around RA,DEC = 0,0\n";
  Stomp::AngularCoordinate fixed_ang(theta, phi,
				     Stomp::AngularCoordinate::Equatorial);
  Stomp::AngularCoordinate start_ang(theta, phi+radius,
				     Stomp::AngularCoordinate::Equatorial);

  // First, we'll try rotating start_ang about fixed_ang by 90 degrees.
  Stomp::AngularCoordinate test_ang;
  start_ang.Rotate(fixed_ang, 90.0, test_ang);
  std::cout << "\t(" << start_ang.RA() << ", " << start_ang.DEC() <<
    ") rotated 90 degrees to (" << test_ang.RA() << ", " <<
    test_ang.DEC() << ")\n";

  // Now, we'll rotate test_ang another 90 degrees.
  test_ang.Rotate(fixed_ang, 90.0);
  std::cout << "\tRotated another 90 degrees to (" << test_ang.RA() <<
    ", " << test_ang.DEC() << ")\n";

  // And finally another 180 degrees to take us back to our starting point.
  test_ang.Rotate(fixed_ang, 180.0);
  std::cout << "\tRotated another 180 degrees to (" << test_ang.RA() <<
    ", " << test_ang.DEC() << ") (should be " << start_ang.RA() << "," <<
    start_ang.DEC() << ")\n";

  // Now do the same tests in Survey coordinates.
  std::cout << "Now Survey coordinates.\n" <<
    "\tRotating around Lambda,Eta = 0,0\n";
  fixed_ang.SetSurveyCoordinates(theta, phi);
  start_ang.SetSurveyCoordinates(theta+radius, phi);

  start_ang.Rotate(fixed_ang, 90.0, test_ang);
  std::cout << "\t(" << start_ang.Lambda() << ", " << start_ang.Eta() <<
    ") rotated 90 degrees to (" << test_ang.Lambda() << ", " <<
    test_ang.Eta() << ")\n";

  // Now, we'll rotate test_ang another 90 degrees.
  test_ang.Rotate(fixed_ang, 90.0);
  std::cout << "\tRotated another 90 degrees to (" << test_ang.Lambda() <<
    ", " << test_ang.Eta() << ")\n";

  // And finally another 180 degrees to take us back to our starting point.
  test_ang.Rotate(fixed_ang, 180.0);
  std::cout << "\tRotated another 180 degrees to (" << test_ang.Lambda() <<
    ", " << test_ang.Eta() << ") (should be " << start_ang.Lambda() << "," <<
    start_ang.Eta() << ")\n";
}

void WeightedAngularCoordinateBasicTests() {
  // Ok, some basic routines to test the extensions to the AngularCoordinate
  // available in WeightedAngularCoordinate.
  std::cout << "\n";
  std::cout << "*********************************************\n";
  std::cout << "*** WeightedAngularCoordinate Basic Tests ***\n";
  std::cout << "*********************************************\n";
  double lambda, eta;

  lambda = 10.0;
  eta = 10.0;

  Stomp::WeightedAngularCoordinate ang(lambda, eta, 1.0,
				       Stomp::AngularCoordinate::Survey);
  ang.SetField("one", 1.0);
  ang.SetField("two", 2.0);
  ang.SetField("three", 3.0);
  ang.SetField("four", 4.0);

  std::cout << "Weight = " << ang.Weight() << "\n";
  std::cout << "Field values:\n";
  std::cout << "\tOne = " << ang.Field("one") << " (1)\n";
  std::cout << "\tTwo = " << ang.Field("two") << " (2)\n";
  std::cout << "\tThree = " << ang.Field("three") << " (3)\n";
  std::cout << "\tFour = " << ang.Field("four") << " (4)\n";
  std::cout << "\tFive = " << ang.Field("five") << " (0)\n";

  Stomp::WeightedAngularCoordinate tmp_ang = ang;
  std::cout << "Copied object:\tWeight = " << tmp_ang.Weight() << "\n";
  std::cout << "\tField values:\n";
  std::cout << "\t\tOne = " << tmp_ang.Field("one") << " (1)\n";
  std::cout << "\t\tTwo = " << tmp_ang.Field("two") << " (2)\n";
  std::cout << "\t\tThree = " << tmp_ang.Field("three") << " (3)\n";
  std::cout << "\t\tFour = " << tmp_ang.Field("four") << " (4)\n";
  std::cout << "\t\tFive = " << tmp_ang.Field("five") << " (0)\n";

  std::cout << "Field copy tests:\n";
  ang.CopyFieldToWeight("two");
  std::cout << "\tAfter moved Field('two') to Weight, Weight = " <<
    ang.Weight() << "\n";
  ang.RestoreOriginalWeight();
  std::cout << "\tAfter restoring Weight, Weight = " << ang.Weight() << "\n";
}

// Define our command line flags
DEFINE_bool(all_angular_coordinate_tests, false, "Run all class unit tests.");
DEFINE_bool(angular_coordinate_basic_tests, false,
            "Run AngularCoordinate basic tests");
DEFINE_bool(angular_coordinate_position_angle_tests, false,
            "Run AngularCoordinate position angle tests");
DEFINE_bool(angular_coordinate_rotation_tests, false,
            "Run AngularCoordinate rotation tests");
DEFINE_bool(weighted_angular_coordinate_basic_tests, false,
            "Run WeightedAngularCoordinate basic tests");

void AngularCoordinateUnitTests(bool run_all_tests) {
  void AngularCoordinateBasicTests();
  void AngularCoordinatePositionAngleTests();
  void AngularCoordinateRotationTests();
  void WeightedAngularCoordinateBasicTests();

  if (run_all_tests) FLAGS_all_angular_coordinate_tests = true;

  // Check that AngularCoordinate transforms from one coordinate system to
  // another work properly and that we can recover the same pixel indices from
  // transformed coordinates.
  if (FLAGS_all_angular_coordinate_tests ||
      FLAGS_angular_coordinate_basic_tests)
    AngularCoordinateBasicTests();

  // Check that the AngularCoordinate position angles are being calculated
  // correctly.
  if (FLAGS_all_angular_coordinate_tests ||
      FLAGS_angular_coordinate_position_angle_tests)
    AngularCoordinatePositionAngleTests();

  // Check that the AngularCoordinate rotations are being calculated
  // correctly.
  if (FLAGS_all_angular_coordinate_tests ||
      FLAGS_angular_coordinate_rotation_tests)
    AngularCoordinateRotationTests();

  // Check WeightedAngularCoordinate extensions to the basic AngularCoordinate.
  if (FLAGS_all_angular_coordinate_tests ||
      FLAGS_weighted_angular_coordinate_basic_tests)
    WeightedAngularCoordinateBasicTests();
}
