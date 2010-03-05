#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>
#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "stomp_map.h"
#include "stomp_footprint.h"

void CircleBoundTests() {
  // Testing the basic pixelization routines.
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

  Stomp::Pixel tmp_pix(ang, Stomp::HPixResolution);
  if (circ.CheckPoint(ang)) {
    std::cout << "Good: CheckPoint inside the bound returns as true (" <<
      Stomp::Footprint::ScorePixel(circ, tmp_pix) << ")\n";
  } else {
    std::cout << "Bad: CheckPoint inside the bound returns as false (" <<
      Stomp::Footprint::ScorePixel(circ, tmp_pix) << ")\n";
  }

  ang.SetSurveyCoordinates(ang.Lambda()+2.0*radius, ang.Eta());
  tmp_pix.SetPixnumFromAng(ang);
  if (circ.CheckPoint(ang)) {
    std::cout << "Bad: CheckPoint outside the bound returns as true (" <<
      Stomp::Footprint::ScorePixel(circ, tmp_pix) << ")\n";
  } else {
    std::cout << "Good: CheckPoint outside the bound returns as false (" <<
      Stomp::Footprint::ScorePixel(circ, tmp_pix) << ")\n";
  }
}

void WedgeBoundTests() {
  // Testing the WedgeBound class
  std::cout << "\n";
  std::cout << "************************\n";
  std::cout << "*** WedgeBound Tests ***\n";
  std::cout << "************************\n";
  Stomp::AngularCoordinate center_ang(20.0, 0.0,
				      Stomp::AngularCoordinate::Survey);
  double radius = 1.0;
  double position_angle_min = 0.0;
  double position_angle_max = 90.0;

  // Start with a wedge that's a quarter of the full circle.
  Stomp::WedgeBound wedge(center_ang, radius,
			  position_angle_min, position_angle_max,
			  Stomp::AngularCoordinate::Survey);
  std::cout << position_angle_min << " - " << position_angle_max <<
    " degrees Wedge Area: " << wedge.Area() << "\n";

  std::cout << "\tBounds:\n" << "\t\tLambda: " << wedge.LambdaMin() << " - " <<
    wedge.LambdaMax() << ", Eta: " << wedge.EtaMin() << " - " <<
    wedge.EtaMax() << "\n";

  std::cout << "\tChecking test points:\n";
  Stomp::AngularCoordinate tmp_ang(center_ang.Lambda()+0.5*radius,
				   center_ang.Eta(),
				   Stomp::AngularCoordinate::Survey);
  tmp_ang.Rotate(center_ang, 0.5*(position_angle_min+position_angle_max));
  Stomp::Pixel tmp_pix(tmp_ang, Stomp::MaxPixelResolution);
  if (wedge.CheckPoint(tmp_ang)) {
    std::cout << "\t\tGood: CheckPoint inside the bound returns as true (" <<
      Stomp::Footprint::ScorePixel(wedge, tmp_pix) << ")\n";
  } else {
    std::cout << "\t\tBad: CheckPoint inside the bound returns as false (" <<
      Stomp::Footprint::ScorePixel(wedge, tmp_pix) << ")\n";
  }

  tmp_ang.Rotate(center_ang, 90.0);
  tmp_pix.SetPixnumFromAng(tmp_ang);
  if (wedge.CheckPoint(tmp_ang)) {
    std::cout << "\t\tBad: CheckPoint outside the bound returns as true (" <<
      Stomp::Footprint::ScorePixel(wedge, tmp_pix) << ")\n";
  } else {
    std::cout << "\t\tGood: CheckPoint outside the bound returns as false (" <<
      Stomp::Footprint::ScorePixel(wedge, tmp_pix) << ")\n";
  }

  // Now we iterate through the other three quadrants to check that we're
  // getting the bounding box right.
  for (int i=0;i<3;i++) {
    position_angle_min += 90.0;
    position_angle_max += 90.0;
    Stomp::WedgeBound new_wedge(center_ang, radius,
				position_angle_min, position_angle_max,
				Stomp::AngularCoordinate::Survey);
    std::cout << "\n" << position_angle_min << " - " << position_angle_max <<
      " degrees Wedge Area: " << wedge.Area() << "\n";

    std::cout << "Bounds:\n" << "\tLambda: " <<
      new_wedge.LambdaMin() << " - " << new_wedge.LambdaMax() << ", Eta: " <<
      new_wedge.EtaMin() << " - " << new_wedge.EtaMax() << "\n";

    std::cout << "\tChecking test points:\n";
    tmp_ang.SetSurveyCoordinates(center_ang.Lambda()+0.5*radius,
				 center_ang.Eta());
    tmp_ang.Rotate(center_ang, 0.5*(position_angle_min+position_angle_max));
    tmp_pix.SetPixnumFromAng(tmp_ang);
    if (new_wedge.CheckPoint(tmp_ang)) {
      std::cout << "\t\tGood: CheckPoint inside the bound returns as true (" <<
	Stomp::Footprint::ScorePixel(new_wedge, tmp_pix) << ")\n";
    } else {
      std::cout << "\t\tBad: CheckPoint inside the bound returns as false (" <<
	Stomp::Footprint::ScorePixel(new_wedge, tmp_pix) << ")\n";
    }

    tmp_ang.Rotate(center_ang, 90.0);
    tmp_pix.SetPixnumFromAng(tmp_ang);
    if (new_wedge.CheckPoint(tmp_ang)) {
      std::cout << "\t\tBad: CheckPoint outside the bound returns as true (" <<
	Stomp::Footprint::ScorePixel(new_wedge, tmp_pix) << ")\n";
    } else {
      std::cout << "\t\tGood: CheckPoint outside the bound returns as false" <<
	" (" << Stomp::Footprint::ScorePixel(new_wedge, tmp_pix) << ")\n";
    }
  }

}

void PolygonBoundTests() {
  // Now we try the polygon-based footprint.
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
  Stomp::Pixel tmp_pix(ang, Stomp::HPixResolution);
  if (poly.CheckPoint(ang)) {
    std::cout << "Good: CheckPoint inside the bound returns as true (" <<
      Stomp::Footprint::ScorePixel(poly, tmp_pix) << ")\n";
  } else {
    std::cout << "Bad: CheckPoint inside the bound returns as false (" <<
      Stomp::Footprint::ScorePixel(poly, tmp_pix) << ")\n";
  }

  ang.SetEquatorialCoordinates(24.0, 0.0);
  tmp_pix.SetPixnumFromAng(ang);
  if (poly.CheckPoint(ang)) {
    std::cout << "Bad: CheckPoint outside the bound returns as true (" <<
      Stomp::Footprint::ScorePixel(poly, tmp_pix) << ")\n";
  } else {
    std::cout << "Good: CheckPoint outside the bound returns as false (" <<
      Stomp::Footprint::ScorePixel(poly, tmp_pix) << ")\n";
  }
}

void FootprintPixelizationTests() {
  // Now we try pixelizing within the Footprint class
  std::cout << "\n";
  std::cout << "***********************\n";
  std::cout << "*** Footprint Tests ***\n";
  std::cout << "***********************\n";

  Stomp::AngularCoordinate ang(20.0, 0.0, Stomp::AngularCoordinate::Survey);
  double radius = 1.0;
  Stomp::CircleBound circ(ang, radius);

  std::cout << "Pixelizing...\n";
  Stomp::Footprint* footprint =
    new Stomp::Footprint(circ, 1.0, Stomp::MaxPixelResolution, false);

  std::cout << "\tStarting resolution level: " <<
    static_cast<int>(footprint->_FindStartingResolutionLevel(circ.Area())) <<
    "\n";

  if (footprint->Pixelize()) {
    std::cout << "Pixelization success!\nArea Check: Real: " <<
      footprint->Area() << ", Pixelized: " <<
      footprint->PixelizedArea() << "\n";
  } else {
    std::cout << "Pixelization failed...\n";
  }

  // Trying again with coarser maximum resolution.
  footprint->SetMaxResolution(2048);
  if (footprint->Pixelize()) {
    std::cout << "Pixelization success!\nArea Check: Real: " <<
      footprint->Area() << ", Pixelized: " <<
      footprint->PixelizedArea() << "\n";
  } else {
    std::cout << "Pixelization failed...\n";
  }
}

void StaticFootprintTests() {
  // Now we try pixelizing within the Footprint class
  std::cout << "\n";
  std::cout << "******************************\n";
  std::cout << "*** Static Footprint Tests ***\n";
  std::cout << "******************************\n";

  Stomp::AngularCoordinate ang(20.0, 0.0, Stomp::AngularCoordinate::Survey);
  double radius = 1.0;
  Stomp::CircleBound circ(ang, radius);
  Stomp::Map stomp_map;

  std::cout << "Pixelizing with static method...\n";
  if (Stomp::Footprint::PixelizeBound(circ, stomp_map, 1.0,
				      Stomp::MaxPixelResolution)) {
    std::cout << "Pixelization success!\nArea Check: Real: " <<
      circ.Area() << ", Pixelized: " << stomp_map.Area() << " (" <<
      stomp_map.Size() << " pixels)\n";
  } else {
    std::cout << "Pixelization failed...\n";
  }

  // Trying again with coarser maximum resolution.
  std::cout << "\nRe-pixelizing at a coarser maximum resolution...\n";
  if (Stomp::Footprint::PixelizeBound(circ, stomp_map, 1.0, 2048)) {
    std::cout << "Pixelization success!\nArea Check: Real: " <<
      circ.Area() << ", Pixelized: " << stomp_map.Area() << " (" <<
      stomp_map.Size() << " pixels)\n";
  } else {
    std::cout << "Pixelization failed...\n";
  }
}


// Define our command line flags
DEFINE_bool(all_footprint_tests, false, "Run all class unit tests.");
DEFINE_bool(circle_bound_tests, false, "Run CircleBound tests");
DEFINE_bool(wedge_bound_tests, false, "Run WedgeBound tests");
DEFINE_bool(polygon_bound_tests, false, "Run PolygonBound tests");
DEFINE_bool(footprint_pixelization_tests, false,
	    "Run Footprint pixelization tests");
DEFINE_bool(static_footprint_tests, false, "Run static Footprint tests");

void FootprintUnitTests(bool run_all_tests) {
  void CircleBoundTests();
  void WedgeBoundTests();
  void PolygonBoundTests();
  void FootprintPixelizationTests();
  void StaticFootprintTests();

  if (run_all_tests) FLAGS_all_footprint_tests = true;

  // Check the CircleBound class.
  if (FLAGS_all_footprint_tests || FLAGS_circle_bound_tests)
    CircleBoundTests();

  // Check the WedgeBound class.
  if (FLAGS_all_footprint_tests || FLAGS_wedge_bound_tests)
    WedgeBoundTests();

  // Check the PolygonBound class.
  if (FLAGS_all_footprint_tests || FLAGS_polygon_bound_tests)
    PolygonBoundTests();

  // Check the Footprint class.
  if (FLAGS_all_footprint_tests || FLAGS_footprint_pixelization_tests)
    FootprintPixelizationTests();

  // Check the static pixelization method in the Footprint class.
  if (FLAGS_all_footprint_tests || FLAGS_static_footprint_tests)
    StaticFootprintTests();
}
