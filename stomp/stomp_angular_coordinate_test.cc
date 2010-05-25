#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <fstream>
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
}

void AngularCoordinateDistanceTests() {
  // Test our routines for calculating the angular distance between points along
  // with those for generating dot products, cross products and great circles.
  std::cout << "\n";
  std::cout << "****************************************\n";
  std::cout << "*** AngularCoordinate Distance tests ***\n";
  std::cout << "****************************************\n";

  Stomp::AngularCoordinate::Sphere sphere = Stomp::AngularCoordinate::Survey;
  Stomp::AngularCoordinate ang_a(0.0, 0.0, sphere);
  Stomp::AngularCoordinate ang_b(45.0, 0.0, sphere);
  Stomp::AngularCoordinate ang_c(60.0, 0.01, sphere);
  Stomp::AngularCoordinate ang_d(0.0, 60.0, sphere);
  std::cout << "\tPoint A: " << ang_a.Lambda() << ", " << ang_a.Eta() << "\n";
  std::cout << "\tPoint B: " << ang_b.Lambda() << ", " << ang_b.Eta() << "\n";
  std::cout << "\tPoint C: " << ang_c.Lambda() << ", " << ang_c.Eta() << "\n";
  std::cout << "\tPoint D: " << ang_d.Lambda() << ", " << ang_d.Eta() << "\n";
  std::cout << "\tA * B: " << ang_a.DotProduct(ang_b) << "\t\t\t B * A: " <<
    ang_b.DotProduct(ang_a) << "\n";
  Stomp::AngularCoordinate ang_n;
  ang_a.GreatCircle(ang_b, ang_n, sphere);
  Stomp::AngularCoordinate ang_m;
  ang_b.GreatCircle(ang_a, ang_m, sphere);

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

void AngularCoordinateUnitSphereTests() {
  // Test that our routines for going back and forth from the Cartesian
  // coordinates are working properly.
  std::cout << "\n";
  std::cout << "*******************************************\n";
  std::cout << "*** AngularCoordinate Unit Sphere tests ***\n";
  std::cout << "*******************************************\n";

  double theta = 60.0, phi = 30.0;

  // Start in Survey coordinates.
  Stomp::AngularCoordinate::Sphere sphere = Stomp::AngularCoordinate::Survey;
  Stomp::AngularCoordinate ang(phi, theta, sphere);
  std::cout << "(Lambda, Eta) = " << ang.Lambda() << "," << ang.Eta() <<
    " -> (" << ang.UnitSphereX() << "," << ang.UnitSphereY() << "," <<
    ang.UnitSphereZ() << ")\n";
  std::cout << "\tSurvey (x,y,z): (" <<
    ang.UnitSphereX(Stomp::AngularCoordinate::Survey) << "," <<
    ang.UnitSphereY(Stomp::AngularCoordinate::Survey) << "," <<
    ang.UnitSphereZ(Stomp::AngularCoordinate::Survey) <<  ")\n";
  std::cout << "\tEquatorial (x,y,z): (" <<
    ang.UnitSphereX(Stomp::AngularCoordinate::Equatorial) << "," <<
    ang.UnitSphereY(Stomp::AngularCoordinate::Equatorial) << "," <<
    ang.UnitSphereZ(Stomp::AngularCoordinate::Equatorial) << ")\n";
  std::cout << "\tGalactic (x,y,z): (" <<
    ang.UnitSphereX(Stomp::AngularCoordinate::Galactic) << "," <<
    ang.UnitSphereY(Stomp::AngularCoordinate::Galactic) << "," <<
    ang.UnitSphereZ(Stomp::AngularCoordinate::Galactic) << ")\n";

  ang.SetEquatorialCoordinates(theta, phi);
  std::cout << "\n(RA, DEC) = " << ang.RA() << "," << ang.DEC() <<
    " -> (" << ang.UnitSphereX() << "," << ang.UnitSphereY() << "," <<
    ang.UnitSphereZ() << ")\n";
  std::cout << "\tSurvey (x,y,z): (" <<
    ang.UnitSphereX(Stomp::AngularCoordinate::Survey) << "," <<
    ang.UnitSphereY(Stomp::AngularCoordinate::Survey) << "," <<
    ang.UnitSphereZ(Stomp::AngularCoordinate::Survey) <<  ")\n";
  std::cout << "\tEquatorial (x,y,z): (" <<
    ang.UnitSphereX(Stomp::AngularCoordinate::Equatorial) << "," <<
    ang.UnitSphereY(Stomp::AngularCoordinate::Equatorial) << "," <<
    ang.UnitSphereZ(Stomp::AngularCoordinate::Equatorial) << ")\n";
  std::cout << "\tGalactic (x,y,z): (" <<
    ang.UnitSphereX(Stomp::AngularCoordinate::Galactic) << "," <<
    ang.UnitSphereY(Stomp::AngularCoordinate::Galactic) << "," <<
    ang.UnitSphereZ(Stomp::AngularCoordinate::Galactic) << ")\n";

  ang.SetGalacticCoordinates(theta, phi);
  std::cout << "\n(GalLon, GalLat) = " << ang.GalLon() << "," << ang.GalLat() <<
    " -> (" << ang.UnitSphereX() << "," << ang.UnitSphereY() << "," <<
    ang.UnitSphereZ() << ")\n";
  std::cout << "\tSurvey (x,y,z): (" <<
    ang.UnitSphereX(Stomp::AngularCoordinate::Survey) << "," <<
    ang.UnitSphereY(Stomp::AngularCoordinate::Survey) << "," <<
    ang.UnitSphereZ(Stomp::AngularCoordinate::Survey) <<  ")\n";
  std::cout << "\tEquatorial (x,y,z): (" <<
    ang.UnitSphereX(Stomp::AngularCoordinate::Equatorial) << "," <<
    ang.UnitSphereY(Stomp::AngularCoordinate::Equatorial) << "," <<
    ang.UnitSphereZ(Stomp::AngularCoordinate::Equatorial) << ")\n";
  std::cout << "\tGalactic (x,y,z): (" <<
    ang.UnitSphereX(Stomp::AngularCoordinate::Galactic) << "," <<
    ang.UnitSphereY(Stomp::AngularCoordinate::Galactic) << "," <<
    ang.UnitSphereZ(Stomp::AngularCoordinate::Galactic) << ")\n";

  std::cout << "\nSetting coordinates based on Cartesian coordinates...\n";
  sphere = Stomp::AngularCoordinate::Survey;
  ang.SetSurveyCoordinates(phi, theta);
  ang.SetUnitSphereCoordinates(ang.UnitSphereX(sphere), ang.UnitSphereY(sphere),
			       ang.UnitSphereZ(sphere), sphere);

  std::cout << "\t(Lambda, Eta) = " << ang.Lambda() << "," << ang.Eta() <<
    " (should be (" << phi << "," << theta << "))\n";
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
  double fixed_theta = 0.0;
  double fixed_phi = 0.0;
  double starting_theta = 0.0;
  double starting_phi = 10.0;
  Stomp::AngularCoordinate::Sphere sphere =
    Stomp::AngularCoordinate::Equatorial;

  Stomp::AngularCoordinate fixed_ang(fixed_theta, fixed_phi, sphere);
  Stomp::AngularCoordinate start_ang(starting_theta, starting_phi, sphere);

  std::cout << "Starting in Equatorial coordinates.\n" <<
    "\tRotating around RA,DEC = " << fixed_ang.RA() << "," <<
    fixed_ang.DEC() << " (" <<
    fixed_ang.UnitSphereX(sphere) << "," <<
    fixed_ang.UnitSphereY(sphere) << "," <<
    fixed_ang.UnitSphereZ(sphere) << ")\n";
  // First, we'll try rotating start_ang about fixed_ang by 90 degrees.
  Stomp::AngularCoordinate test_ang;
  start_ang.Rotate(fixed_ang, 90.0, test_ang, sphere);
  std::cout << "\t(" << start_ang.RA() << ", " << start_ang.DEC() <<
    ") rotated 90 degrees to (" << test_ang.RA() << ", " <<
    test_ang.DEC() << ")\n\t\t(" <<
    test_ang.UnitSphereX(sphere) << "," <<
    test_ang.UnitSphereY(sphere) << "," <<
    test_ang.UnitSphereZ(sphere) << ")\n";

  // Now, we'll rotate test_ang another 90 degrees.
  test_ang.Rotate(fixed_ang, 90.0, sphere);
  std::cout << "\tRotated another 90 degrees to (" << test_ang.RA() <<
    ", " << test_ang.DEC() << ")\n\t\t(" <<
    test_ang.UnitSphereX(sphere) << "," <<
    test_ang.UnitSphereY(sphere) << "," <<
    test_ang.UnitSphereZ(sphere) << ")\n";

  // And finally another 180 degrees to take us back to our starting point.
  test_ang.Rotate(fixed_ang, 180.0, sphere);
  std::cout << "\tRotated another 180 degrees to (" << test_ang.RA() <<
    ", " << test_ang.DEC() << ") (should be " << start_ang.RA() << "," <<
    start_ang.DEC() << ")\n\t\t(" <<
    test_ang.UnitSphereX(sphere) << "," <<
    test_ang.UnitSphereY(sphere) << "," <<
    test_ang.UnitSphereZ(sphere) << ")\n";


  // Now do the same tests in Survey coordinates.
  sphere = Stomp::AngularCoordinate::Survey;
  fixed_ang.SetSurveyCoordinates(fixed_phi, fixed_theta);
  start_ang.SetSurveyCoordinates(starting_phi, starting_theta);
  std::cout << "\nNow Survey coordinates.\n" <<
    "\tRotating around Lambda,Eta = " << fixed_ang.Lambda() << "," <<
    fixed_ang.Eta() << " (" <<
    fixed_ang.UnitSphereX(sphere) << "," <<
    fixed_ang.UnitSphereY(sphere) << "," <<
    fixed_ang.UnitSphereZ(sphere) << ")\n";

  start_ang.Rotate(fixed_ang, 90.0, test_ang, sphere);
  std::cout << "\t(" << start_ang.Lambda() << ", " << start_ang.Eta() <<
    ") rotated 90 degrees to (" << test_ang.Lambda() << ", " <<
    test_ang.Eta() << ")\n\t\t(" <<
    test_ang.UnitSphereX(sphere) << "," <<
    test_ang.UnitSphereY(sphere) << "," <<
    test_ang.UnitSphereZ(sphere) << ")\n";

  // Now, we'll rotate test_ang another 90 degrees.
  test_ang.Rotate(fixed_ang, 90.0, sphere);
  std::cout << "\tRotated another 90 degrees to (" << test_ang.Lambda() <<
    ", " << test_ang.Eta() << ")\n\t\t(" <<
    test_ang.UnitSphereX(sphere) << "," <<
    test_ang.UnitSphereY(sphere) << "," <<
    test_ang.UnitSphereZ(sphere) << ")\n";

  // And finally another 180 degrees to take us back to our starting point.
  test_ang.Rotate(fixed_ang, 180.0, sphere);
  std::cout << "\tRotated another 180 degrees to (" << test_ang.Lambda() <<
    ", " << test_ang.Eta() << ") (should be " << start_ang.Lambda() << "," <<
    start_ang.Eta() << ")\n\t\t(" <<
    test_ang.UnitSphereX(sphere) << "," <<
    test_ang.UnitSphereY(sphere) << "," <<
    test_ang.UnitSphereZ(sphere) << ")\n";
}

void AngularVectorIOTests() {
  // Check to make sure that our static methods for converting between
  // std::vectors and AngularVectors are working properly..
  std::cout << "\n";
  std::cout << "*******************************\n";
  std::cout << "*** AngularVector I/O Check ***\n";
  std::cout << "*******************************\n\n";

  // Start by making some RA & DEC vectors.
  std::vector<double> raVec, decVec;

  raVec.push_back(0.0);
  raVec.push_back(20.0);
  raVec.push_back(90.0);

  decVec.push_back(0.0);
  decVec.push_back(50.0);
  decVec.push_back(-90.0);

  // Now let's make an AngularVector out of those coordinates.
  Stomp::AngularCoordinate::Sphere sphere =
    Stomp::AngularCoordinate::Equatorial;
  Stomp::AngularVector ang;

  if (Stomp::AngularCoordinate::ToAngularVector(raVec, decVec, ang, sphere)) {
    std::cout << "\tConverted doubles vectors to AngularVector...\n\n";
  } else {
    std::cout << "\tFailed to convert doubles vectors to AngularVector...\n\n";
  }

  // And generate new doubles vectors from that AngularVector
  std::vector<double> raVec_out, decVec_out;
  if (Stomp::AngularCoordinate::FromAngularVector(ang, raVec_out,
						  decVec_out, sphere)) {
    std::cout << "\tConverted AngularVector to doubles vectors...\n\n";
  } else {
    std::cout << "\tFailed to convert AngularVector to doubles vectors...\n\n";
  }

  // And finish up with some visual output to verify that things are working.
  std::cout << "Coordinate check:\n" <<
    "\tInput vector\tAngularVector\tOutputVector\n" <<
    "\t------------\t-------------\t------------\n";
  for (uint8_t i=0;i<raVec.size();i++) {
    std::cout << "\t" << raVec[i] << "," << decVec[i] << "\t\t" <<
      ang[i].RA() << "," << ang[i].DEC() << "\t\t" <<
      raVec_out[i] << "," << decVec_out[i] << "\n";
  }
}

void AngularVectorFileIOTests() {
  // Check to make sure that our static methods for reading ASCII files into
  // AngularVectors are working properly..
  std::cout << "\n";
  std::cout << "************************************\n";
  std::cout << "*** AngularVector File I/O Check ***\n";
  std::cout << "************************************\n\n";

  // Start by making some RA & DEC vectors.
  std::vector<double> raVec, decVec;

  raVec.push_back(0.0);
  raVec.push_back(20.0);
  raVec.push_back(90.0);

  decVec.push_back(0.0);
  decVec.push_back(50.0);
  decVec.push_back(-90.0);

  // Now we'll write them to the simplest file to parse
  std::string test_file = "AngularVectorTest.dat";
  std::ofstream test_file_str(test_file.c_str());

  if (test_file_str.is_open()) {
    for (uint8_t i=0;i<raVec.size();i++) {
      test_file_str << raVec[i] << " " << decVec[i] << "\n";
    }
  }

  test_file_str.close();

  // Now we attempt to read in those coordinates.
  Stomp::AngularCoordinate::Sphere sphere =
    Stomp::AngularCoordinate::Equatorial;
  Stomp::AngularVector ang;

  if (Stomp::AngularCoordinate::ToAngularVector(test_file, ang, sphere)) {
    std::cout << "\tRead " << test_file << " to AngularVector...\n\n";
  } else {
    std::cout << "\tFailed to read " << test_file << " to AngularVector...\n\n";
  }

  // Check to make sure that our coordinates were read in correctly.
  std::cout << "Coordinate check:\n" <<
    "\tFile coordinate\tAngularVector\n" <<
    "\t---------------\t-------------\n";
  for (uint8_t i=0;i<raVec.size();i++) {
    std::cout << "\t" << raVec[i] << "," << decVec[i] << "\t\t" <<
      ang[i].RA() << "," << ang[i].DEC() << "\n";
  }
  std::cout << "\n";

  // Now a trickier version where we've mixed up the coordinates and added
  // some intervening columns.
  std::ofstream test_file2_str(test_file.c_str());

  if (test_file2_str.is_open()) {
    for (uint8_t i=0;i<raVec.size();i++) {
      test_file2_str << "0.0 0.0 " << decVec[i] << " 0.0 0.0 0.0 " <<
	raVec[i] << " 0.0 0.0 0.0 0.0\n";
    }
  }

  test_file2_str.close();

  if (Stomp::AngularCoordinate::ToAngularVector(test_file, ang, sphere, 6 ,2)) {
    std::cout << "\tRead " << test_file << " to AngularVector...\n\n";
  } else {
    std::cout << "\tFailed to read " << test_file << " to AngularVector...\n\n";
  }

  // Check to make sure that our coordinates were read in correctly.
  std::cout << "Coordinate check:\n" <<
    "\tFile coordinate\tAngularVector\n" <<
    "\t---------------\t-------------\n";
  for (uint8_t i=0;i<raVec.size();i++) {
    std::cout << "\t" << raVec[i] << "," << decVec[i] << "\t\t" <<
      ang[i].RA() << "," << ang[i].DEC() << "\n";
  }
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

void WAngularVectorFileIOTests() {
  // Check to make sure that our static methods for reading ASCII files into
  // AngularVectors are working properly..
  std::cout << "\n";
  std::cout << "*************************************\n";
  std::cout << "*** WAngularVector File I/O Check ***\n";
  std::cout << "*************************************\n\n";

  // Start by making some RA & DEC vectors.
  std::vector<double> raVec, decVec, weightVec;

  raVec.push_back(0.0);
  raVec.push_back(20.0);
  raVec.push_back(90.0);

  decVec.push_back(0.0);
  decVec.push_back(50.0);
  decVec.push_back(-90.0);

  weightVec.push_back(1.0);
  weightVec.push_back(2.0);
  weightVec.push_back(3.0);

  // Now we'll write them to the simplest file to parse
  std::string test_file = "WAngularVectorTest.dat";
  std::ofstream test_file_str(test_file.c_str());

  if (test_file_str.is_open()) {
    for (uint8_t i=0;i<raVec.size();i++) {
      test_file_str << raVec[i] << " " << decVec[i] <<
	" " << weightVec[i] << "\n";
    }
  }

  test_file_str.close();

  // Now we attempt to read in those coordinates.
  Stomp::AngularCoordinate::Sphere sphere =
    Stomp::AngularCoordinate::Equatorial;
  Stomp::WAngularVector w_ang;

  if (Stomp::WeightedAngularCoordinate::ToWAngularVector(test_file, w_ang,
							 sphere, false,
							 0, 1, 2)) {
    std::cout << "Read " << test_file << " to WAngularVector...\n";
  } else {
    std::cout << "Failed to read " << test_file << " to WAngularVector...\n";
  }

  // Check to make sure that our coordinates were read in correctly.
  std::cout << "Coordinate check:\n" <<
    "\tFile coordinate\tWAngularVector\n" <<
    "\t---------------\t--------------\n";
  for (uint8_t i=0;i<raVec.size();i++) {
    std::cout << "\t" <<
      raVec[i] << "," << decVec[i] << "," << weightVec[i] << "\t\t" <<
      w_ang[i].RA() << "," << w_ang[i].DEC() << "," <<
      w_ang[i].Weight() << "\n";
  }
  std::cout << "\n";


  // Now a trickier version where we've mixed up the coordinates and added
  // some intervening columns.
  std::ofstream test_file2_str(test_file.c_str());

  if (test_file2_str.is_open()) {
    for (uint8_t i=0;i<raVec.size();i++) {
      test_file2_str << "0.0 0.0 " << decVec[i] << " 0.0 0.0 0.0 " <<
	raVec[i] << " 0.0 0.0 0.0 " << weightVec[i] << " 0.0 0.0\n";
    }
  }

  test_file2_str.close();

  if (Stomp::WeightedAngularCoordinate::ToWAngularVector(test_file, w_ang,
							 sphere, false,
							 6 , 2, 10)) {
    std::cout << "Read " << test_file << " to WAngularVector...\n";
  } else {
    std::cout << "Failed to read " << test_file << " to WAngularVector...\n";
  }

  // Check to make sure that our coordinates were read in correctly.
  std::cout << "Coordinate check:\n" <<
    "\tFile coordinate\tWAngularVector\n" <<
    "\t---------------\t--------------\n";
  for (uint8_t i=0;i<raVec.size();i++) {
    std::cout << "\t" <<
      raVec[i] << "," << decVec[i] << "," << weightVec[i] << "\t\t" <<
      w_ang[i].RA() << "," << w_ang[i].DEC() << "," <<
      w_ang[i].Weight() << "\n";
  }
  std::cout << "\n";


  // Finally, we add in the final wrinkle, bringing in the Field values.
  std::ofstream test_file3_str(test_file.c_str());

  if (test_file3_str.is_open()) {
    for (uint8_t i=0;i<raVec.size();i++) {
      test_file3_str << "1.0 2.0 " << decVec[i] << " 3.0 4.0 5.0 " <<
	raVec[i] << " 6.0 7.0 8.0 " << weightVec[i] << " 9.0 10.0\n";
    }
  }

  test_file3_str.close();

  Stomp::FieldColumnDict field_names;
  field_names["One"] = 0;
  // field_names["Two"] = 1;
  // field_names["Three"] = 3;
  // field_names["Four"] = 4;
  field_names["Five"] = 5;
  // field_names["Six"] = 7;
  // field_names["Seven"] = 8;
  // field_names["Eight"] = 9;
  field_names["Nine"] = 11;
  // field_names["Ten"] = 12;

  if (Stomp::WeightedAngularCoordinate::ToWAngularVector(test_file, w_ang,
							 field_names, sphere,
							 false, 6 , 2, 10)) {
    std::cout << "Read " << test_file << " to WAngularVector...\n";
  } else {
    std::cout << "Failed to read " << test_file << " to WAngularVector...\n";
  }

  // Check to make sure that our coordinates were read in correctly.
  std::cout << "Coordinate check:\n" <<
    "\tFile coordinate\tWAngularVector\n" <<
    "\t---------------\t--------------\n";
  for (uint8_t i=0;i<raVec.size();i++) {
    std::cout << "\t" <<
      raVec[i] << "," << decVec[i] << "," << weightVec[i] << "\t\t" <<
      w_ang[i].RA() << "," << w_ang[i].DEC() << "," << w_ang[i].Weight() <<
      ",Fields: ";
    for (Stomp::FieldIterator iter=w_ang[i].FieldBegin();
	 iter!=w_ang[i].FieldEnd();++iter) {
      std::cout << iter->first << ":" << w_ang[i].Field(iter->first) << ", ";
    }
    std::cout << "\n";
  }
}

// Define our command line flags
DEFINE_bool(all_angular_coordinate_tests, false, "Run all class unit tests.");
DEFINE_bool(angular_coordinate_basic_tests, false,
            "Run AngularCoordinate basic tests");
DEFINE_bool(angular_coordinate_distance_tests, false,
            "Run AngularCoordinate angular distance tests");
DEFINE_bool(angular_coordinate_unit_sphere_tests, false,
            "Run AngularCoordinate Cartesian coordinate tests");
DEFINE_bool(angular_coordinate_position_angle_tests, false,
            "Run AngularCoordinate position angle tests");
DEFINE_bool(angular_coordinate_rotation_tests, false,
            "Run AngularCoordinate rotation tests");
DEFINE_bool(angular_vector_io_tests, false,
            "Run AngularVector I/O tests ");
DEFINE_bool(angular_vector_file_io_tests, false,
            "Run AngularVector file I/O tests ");
DEFINE_bool(weighted_angular_coordinate_basic_tests, false,
            "Run WeightedAngularCoordinate basic tests");
DEFINE_bool(wangular_vector_file_io_tests, false,
            "Run WAngularVector file I/O tests ");

void AngularCoordinateUnitTests(bool run_all_tests) {
  void AngularCoordinateBasicTests();
  void AngularCoordinateDistanceTests();
  void AngularCoordinateUnitSphereTests();
  void AngularCoordinatePositionAngleTests();
  void AngularCoordinateRotationTests();
  void AngularVectorIOTests();
  void AngularVectorFileIOTests();
  void WeightedAngularCoordinateBasicTests();
  void WAngularVectorFileIOTests();

  if (run_all_tests) FLAGS_all_angular_coordinate_tests = true;

  // Check that AngularCoordinate transforms from one coordinate system to
  // another work properly and that we can recover the same pixel indices from
  // transformed coordinates.
  if (FLAGS_all_angular_coordinate_tests ||
      FLAGS_angular_coordinate_basic_tests)
    AngularCoordinateBasicTests();

  // Check that our routines for calculating angular distances, dot and cross
  // products are working.
  if (FLAGS_all_angular_coordinate_tests ||
      FLAGS_angular_coordinate_distance_tests)
    AngularCoordinateDistanceTests();

  // Check that our routines for calculating and setting the Cartesian
  // coordinates are are working properly
  if (FLAGS_all_angular_coordinate_tests ||
      FLAGS_angular_coordinate_unit_sphere_tests)
    AngularCoordinateUnitSphereTests();

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

  // Check that the conversions between vector<double> and AngularVector are
  // working properly
  if (FLAGS_all_angular_coordinate_tests || FLAGS_angular_vector_io_tests)
    AngularVectorIOTests();

  // Check that the methods to read ASCII files into AngularVectors are
  // working properly
  if (FLAGS_all_angular_coordinate_tests || FLAGS_angular_vector_file_io_tests)
    AngularVectorFileIOTests();

  // Check WeightedAngularCoordinate extensions to the basic AngularCoordinate.
  if (FLAGS_all_angular_coordinate_tests ||
      FLAGS_weighted_angular_coordinate_basic_tests)
    WeightedAngularCoordinateBasicTests();

  // Check that the methods to read ASCII files into WAngularVectors are
  // working properly
  if (FLAGS_all_angular_coordinate_tests ||
      FLAGS_wangular_vector_file_io_tests)
    WAngularVectorFileIOTests();
}
