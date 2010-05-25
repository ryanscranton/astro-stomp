#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>
#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "stomp_geometry.h"

void CircleBoundTests() {
  // Testing the CircleBound class
  std::cout << "\n";
  std::cout << "*************************\n";
  std::cout << "*** CircleBound Tests ***\n";
  std::cout << "*************************\n";
  Stomp::AngularCoordinate ang(20.0, 0.0, Stomp::AngularCoordinate::Survey);
  double radius = 1.0;

  Stomp::CircleBound circ(ang, radius);
  std::cout << "Circle Area: " << circ.Area() << " (" <<
    (1.0 - cos(radius*Stomp::DegToRad))*
    2.0*Stomp::Pi*Stomp::StradToDeg << ")\n";

  std::cout << "Bounds:\n" << "\tLambda: " << circ.LambdaMin() << " - " <<
    circ.LambdaMax() << ", Eta: " << circ.EtaMin() << " - " <<
    circ.EtaMax() << "\n";

  if (circ.CheckPoint(ang)) {
    std::cout << "Good: CheckPoint inside the bound returns as true\n";
  } else {
    std::cout << "Bad: CheckPoint inside the bound returns as false\n";
  }

  ang.SetSurveyCoordinates(ang.Lambda()+2.0*radius, ang.Eta());
  if (circ.CheckPoint(ang)) {
    std::cout << "Bad: CheckPoint outside the bound returns as true\n";
  } else {
    std::cout << "Good: CheckPoint outside the bound returns as false\n";
  }
}

void AnnulusBoundTests() {
  // Testing the AnnulusBound class
  std::cout << "\n";
  std::cout << "**************************\n";
  std::cout << "*** AnnulusBound Tests ***\n";
  std::cout << "**************************\n";
  Stomp::AngularCoordinate ang(20.0, 0.0, Stomp::AngularCoordinate::Survey);
  double min_radius = 0.5;
  double max_radius = 1.0;

  Stomp::AnnulusBound annulus(ang, min_radius, max_radius);
  std::cout << "Annulus Area: " << annulus.Area() << " (" <<
    (cos(min_radius*Stomp::DegToRad) - cos(max_radius*Stomp::DegToRad))*
    2.0*Stomp::Pi*Stomp::StradToDeg << ")\n";

  Stomp::AngularBin angular_bin(min_radius, max_radius);
  Stomp::AnnulusBound alt_annulus(ang, angular_bin);
  std::cout << "Alternate constructor annulus Area: " <<
    annulus.Area() << " (" <<
    (cos(min_radius*Stomp::DegToRad) - cos(max_radius*Stomp::DegToRad))*
    2.0*Stomp::Pi*Stomp::StradToDeg << ")\n";
  
  std::cout << "Bounds:\n" << "\tLambda: " << annulus.LambdaMin() << " - " <<
    annulus.LambdaMax() << ", Eta: " << annulus.EtaMin() << " - " <<
    annulus.EtaMax() << "\n";

  if (annulus.CheckPoint(ang)) {
    std::cout << "Bad: CheckPoint at the center returns as true\n";
  } else {
    std::cout << "Good: CheckPoint at the center returns as false\n";
  }

  ang.SetSurveyCoordinates(ang.Lambda()+0.5*(min_radius+max_radius), ang.Eta());
  if (annulus.CheckPoint(ang)) {
    std::cout << "Good: CheckPoint inside the bound returns as true\n";
  } else {
    std::cout << "Bad: CheckPoint inside the bound returns as false\n";
  }

  ang.SetSurveyCoordinates(ang.Lambda()+2.0*max_radius, ang.Eta());
  if (annulus.CheckPoint(ang)) {
    std::cout << "Bad: CheckPoint outside the bound returns as true\n";
  } else {
    std::cout << "Good: CheckPoint outside the bound returns as false\n";
  }
}

void WedgeBoundTests() {
  // Testing the WedgeBound class
  std::cout << "\n";
  std::cout << "************************\n";
  std::cout << "*** WedgeBound Tests ***\n";
  std::cout << "************************\n";
  // First, we test a wedge in Survey coordinates.
  Stomp::AngularCoordinate::Sphere sphere = Stomp::AngularCoordinate::Survey;
  Stomp::AngularCoordinate center_ang(20.0, 0.0, sphere);
  double radius = 1.0;
  std::cout << "Testing wedges at Lambda,Eta = " << center_ang.Lambda() <<
    "," << center_ang.Eta() << " with " << radius << " degree radius...\n";

  for (int i=0;i<4;i++) {
    double position_angle_min = 0.0 + i*90.0;
    double position_angle_max = 90.0 + i*90.0;
    Stomp::WedgeBound wedge(center_ang, radius,
			    position_angle_min, position_angle_max, sphere);
    std::cout << "\n" << position_angle_min << " - " << position_angle_max <<
      " degrees Wedge Area: " << wedge.Area() << "\n";

    std::cout << "Bounds:\n" << "\tLambda: " <<
      wedge.LambdaMin() << " - " << wedge.LambdaMax() << ", Eta: " <<
      wedge.EtaMin() << " - " << wedge.EtaMax() << "\n";

    std::cout << "Checking test points:\n";
    Stomp::AngularCoordinate tmp_ang;
    tmp_ang.SetSurveyCoordinates(center_ang.Lambda()+0.5*radius,
				 center_ang.Eta());
    tmp_ang.Rotate(center_ang, 0.5*(position_angle_min+position_angle_max),
		   sphere);
    if (wedge.CheckPoint(tmp_ang)) {
      std::cout << "\tGood: CheckPoint inside the bound (Lambda,Eta = " <<
	tmp_ang.Lambda() << "," << tmp_ang.Eta() << "; PA = " <<
	center_ang.PositionAngle(tmp_ang, sphere) << ") returns as true\n";
    } else {
      std::cout << "\tBad: CheckPoint inside the bound (Lambda,Eta = " <<
	tmp_ang.Lambda() << "," << tmp_ang.Eta() << "; PA = " <<
	center_ang.PositionAngle(tmp_ang, sphere) << ") returns as false\n";
    }

    tmp_ang.Rotate(center_ang, 90.0, sphere);
    if (wedge.CheckPoint(tmp_ang)) {
      std::cout << "\tBad: CheckPoint outside the bound (Lambda,Eta = " <<
	tmp_ang.Lambda() << "," << tmp_ang.Eta() << "; PA = " <<
	center_ang.PositionAngle(tmp_ang, sphere) << ") returns as true\n";
    } else {
      std::cout << "\tGood: CheckPoint outside the bound (Lambda,Eta = " <<
	tmp_ang.Lambda() << "," << tmp_ang.Eta() << "; PA = " <<
	center_ang.PositionAngle(tmp_ang, sphere) << ") returns as false\n";
    }
  }

  // Now we test a similar set of wedges in Equatorial coordinates.
  sphere = Stomp::AngularCoordinate::Equatorial;
  center_ang.SetEquatorialCoordinates(0.0, 20.0);
  std::cout << "\n\nTesting wedges at RA,DEC = " << center_ang.RA() <<
    "," << center_ang.DEC() << " with " << radius << " degree radius...\n";

  for (int i=0;i<4;i++) {
    double position_angle_min = 0.0 + i*90.0;
    double position_angle_max = 90.0 + i*90.0;
    Stomp::WedgeBound wedge(center_ang, radius,
			    position_angle_min, position_angle_max, sphere);
    std::cout << "\n" << position_angle_min << " - " << position_angle_max <<
      " degrees Wedge Area: " << wedge.Area() << "\n";

    std::cout << "Bounds:\n" << "\tLambda: " <<
      wedge.LambdaMin() << " - " << wedge.LambdaMax() << ", Eta: " <<
      wedge.EtaMin() << " - " << wedge.EtaMax() << "\n";

    std::cout << "Checking test points:\n";
    Stomp::AngularCoordinate tmp_ang;
    tmp_ang.SetEquatorialCoordinates(center_ang.RA(),
				     center_ang.DEC()+0.5*radius);
    tmp_ang.Rotate(center_ang, 0.5*(position_angle_min+position_angle_max),
		   sphere);
    if (wedge.CheckPoint(tmp_ang)) {
      std::cout << "\tGood: CheckPoint inside the bound (RA,DEC = " <<
	tmp_ang.RA() << "," << tmp_ang.DEC() << "; PA = " <<
	center_ang.PositionAngle(tmp_ang, sphere) << ") returns as true\n";
    } else {
      std::cout << "\tBad: CheckPoint inside the bound (RA,DEC = " <<
	tmp_ang.RA() << "," << tmp_ang.DEC() << "; PA = " <<
	center_ang.PositionAngle(tmp_ang, sphere) << ") returns as false\n";
    }

    tmp_ang.Rotate(center_ang, 90.0, sphere);
    if (wedge.CheckPoint(tmp_ang)) {
      std::cout << "\tBad: CheckPoint outside the bound (RA,DEC = " <<
	tmp_ang.RA() << "," << tmp_ang.DEC() << "; PA = " <<
	center_ang.PositionAngle(tmp_ang, sphere) << ") returns as true\n";
    } else {
      std::cout << "\tGood: CheckPoint outside the bound (RA,DEC = " <<
	tmp_ang.RA() << "," << tmp_ang.DEC() << "; PA = " <<
	center_ang.PositionAngle(tmp_ang, sphere) << ") returns as false\n";
    }
  }

}

void PolygonBoundTests() {
  // Test the PolygonBound class
  std::cout << "\n";
  std::cout << "**************************\n";
  std::cout << "*** PolygonBound Tests ***\n";
  std::cout << "**************************\n";

  Stomp::AngularCoordinate ang;
  Stomp::AngularVector angVec;
  ang.SetEquatorialCoordinates(17.0,3.0);
  angVec.push_back(ang);
  ang.SetEquatorialCoordinates(17.0,-3.0);
  angVec.push_back(ang);
  ang.SetEquatorialCoordinates(23.0,-3.0);
  angVec.push_back(ang);
  ang.SetEquatorialCoordinates(23.0,3.0);
  angVec.push_back(ang);

  Stomp::PolygonBound poly(angVec);
  std::cout << "Polygon Area: " << poly.Area() << "\n";

  std::cout << "Bounds:\n" << "\tLambda: " <<
    poly.LambdaMin() << " - " << poly.LambdaMax() << ", Eta: " <<
    poly.EtaMin() << " - " << poly.EtaMax() << "\n";

  ang.SetEquatorialCoordinates(20.0, 0.0);
  if (poly.CheckPoint(ang)) {
    std::cout << "Good: CheckPoint inside the bound returns as true\n";
  } else {
    std::cout << "Bad: CheckPoint inside the bound returns as false\n";
  }

  ang.SetEquatorialCoordinates(24.0, 0.0);
  if (poly.CheckPoint(ang)) {
    std::cout << "Bad: CheckPoint outside the bound returns as true\n";
  } else {
    std::cout << "Good: CheckPoint outside the bound returns as false\n";
  }
}

// Define our command line flags
DEFINE_bool(all_geometry_tests, false, "Run all class unit tests.");
DEFINE_bool(circle_bound_tests, false, "Run CircleBound tests");
DEFINE_bool(annulus_bound_tests, false, "Run AnnulusBound tests");
DEFINE_bool(wedge_bound_tests, false, "Run WedgeBound tests");
DEFINE_bool(polygon_bound_tests, false, "Run PolygonBound tests");

void GeometryUnitTests(bool run_all_tests) {
  void CircleBoundTests();
  void AnnulusBoundTests();
  void WedgeBoundTests();
  void PolygonBoundTests();

  if (run_all_tests) FLAGS_all_geometry_tests = true;

  // Check the CircleBound class.
  if (FLAGS_all_geometry_tests || FLAGS_circle_bound_tests)
    CircleBoundTests();

  // Check the AnnulusBound class.
  if (FLAGS_all_geometry_tests || FLAGS_annulus_bound_tests)
    AnnulusBoundTests();

  // Check the WedgeBound class.
  if (FLAGS_all_geometry_tests || FLAGS_wedge_bound_tests)
    WedgeBoundTests();

  // Check the PolygonBound class.
  if (FLAGS_all_geometry_tests || FLAGS_polygon_bound_tests)
    PolygonBoundTests();
}
