#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>
#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "stomp_footprint_test.h"

void CircleFootprintTests() {
  // Testing the basic pixelization routines.
  std::cout << "\n";
  std::cout << "********************************\n";
  std::cout << "*** Circular Footprint Tests ***\n";
  std::cout << "********************************\n";
  Stomp::AngularCoordinate ang(20.0, 0.0, Stomp::AngularCoordinate::Equatorial);

  Stomp::CircleBound* circ = new Stomp::CircleBound(ang, 3.0, 1.0);
  std::cout << "Circle Area: " << circ->Area() << " (" <<
    (1.0 - cos(3.0*Stomp::DegToRad))*
    2.0*Stomp::Pi*Stomp::StradToDeg << ")\n";
  std::cout << "Pixelizing...\n";
  std::cout << "\tStarting resolution level: " <<
    static_cast<int>(circ->FindStartingResolutionLevel()) << "\n";

  if (circ->FindXYBounds(circ->FindStartingResolutionLevel())) {
    std::cout << "\t\tFound X-Y bounds: " <<
        circ->XMin() << " - " << circ->XMax() << ", " <<
        circ->YMin() << " - " << circ->YMax() << "\n";
  } else {
    std::cout << "\t\tFindXYBounds failed...\n";
  }

  if (circ->Pixelize()) {
    std::cout << "Pixelization success!\nArea Check: Real: " <<
        circ->Area() << ", Pixelized: " << circ->PixelizedArea() << "\n";
  } else {
    std::cout << "Pixelization failed...\n";
  }

  // Trying again with coarser maximum resolution.
  circ->SetMaxResolution(2048);
  if (circ->Pixelize()) {
    std::cout << "Pixelization success!\nArea Check: Real: " <<
        circ->Area() << ", Pixelized: " << circ->PixelizedArea() << "\n";
  } else {
    std::cout << "Pixelization failed...\n";
  }
}

void WedgeFootprintTests() {
  // Testing the basic pixelization routines.
  std::cout << "\n";
  std::cout << "*****************************\n";
  std::cout << "*** Wedge Footprint Tests ***\n";
  std::cout << "*****************************\n";
  Stomp::AngularCoordinate ang(20.0, 0.0, Stomp::AngularCoordinate::Equatorial);
  double radius = 1.0;
  double position_angle_min = 0.0;
  double position_angle_max = 90.0;

  // Start with a wedge that's a quarter of the full circle.
  Stomp::WedgeBound* wedge =
    new Stomp::WedgeBound(ang, radius, position_angle_min, position_angle_max,
			  1.0, Stomp::AngularCoordinate::Equatorial);
  std::cout << position_angle_min << " - " << position_angle_max <<
    " degrees Wedge Area: " << wedge->Area() << "\n";
  std::cout << "Pixelizing...\n";
  std::cout << "\tStarting resolution level: " <<
    static_cast<int>(wedge->FindStartingResolutionLevel()) << "\n";

  if (wedge->FindXYBounds(wedge->FindStartingResolutionLevel())) {
    std::cout << "\t\tFound X-Y bounds: " <<
        wedge->XMin() << " - " << wedge->XMax() << ", " <<
        wedge->YMin() << " - " << wedge->YMax() << "\n";
  } else {
    std::cout << "\t\tFindXYBounds failed...\n";
  }

  if (wedge->Pixelize()) {
    std::cout << "Pixelization success!\nArea Check: Real: " <<
      wedge->Area() << ", Pixelized: " << wedge->PixelizedArea() << "\n";
  } else {
    std::cout << "Pixelization failed...\n";
  }

  delete wedge;

  // Now we iterate through the other three quadrants to check that we're
  // getting the bounding box right.
  for (int i=0;i<4;i++) {
    position_angle_min += 90.0;
    position_angle_max += 90.0;
    wedge =
      new Stomp::WedgeBound(ang, radius, position_angle_min, position_angle_max,
			    1.0, Stomp::AngularCoordinate::Equatorial);
    std::cout << "\n" << position_angle_min << " - " << position_angle_max <<
      " degrees Wedge Area: " << wedge->Area() << "\n";
    std::cout << "Pixelizing...\n";
    std::cout << "\tStarting resolution level: " <<
      static_cast<int>(wedge->FindStartingResolutionLevel()) << "\n";
    
    if (wedge->FindXYBounds(wedge->FindStartingResolutionLevel())) {
      std::cout << "\t\tFound X-Y bounds: " <<
        wedge->XMin() << " - " << wedge->XMax() << ", " <<
        wedge->YMin() << " - " << wedge->YMax() << "\n";
    } else {
      std::cout << "\t\tFindXYBounds failed...\n";
    }

    if (wedge->Pixelize()) {
      std::cout << "Pixelization success!\nArea Check: Real: " <<
	wedge->Area() << ", Pixelized: " << wedge->PixelizedArea() << "\n";
    } else {
      std::cout << "Pixelization failed...\n";
    }
  }

}

void PolygonFootprintTests() {
  // Now we try the polygon-based footprint.
  std::cout << "\n";
  std::cout << "*******************************\n";
  std::cout << "*** Polygon Footprint Tests ***\n";
  std::cout << "*******************************\n";

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

  Stomp::PolygonBound* poly = new Stomp::PolygonBound(angVec,1.0);
  std::cout << "Polygon Area: " << poly->Area() << "\n";
  std::cout << "Pixelizing...\n";
  std::cout << "\tStarting resolution level: " <<
    static_cast<int>(poly->FindStartingResolutionLevel()) << "\n";
  if (poly->FindXYBounds(poly->FindStartingResolutionLevel())) {
    std::cout << "\t\tFound X-Y bounds: " <<
        poly->XMin() << " - " << poly->XMax() << ", " <<
        poly->YMin() << " - " << poly->YMax() << "\n";
  } else {
    std::cout << "\t\tFindXYBounds failed...\n";
  }
  if (poly->Pixelize()) {
    std::cout << "Pixelization success!\nArea Check: Real: " <<
        poly->Area() << ", Pixelized: " << poly->PixelizedArea() << "\n";
  } else {
    std::cout << "Pixelization failed...\n";
  }

  // Trying again with coarser maximum resolution.
  poly->SetMaxResolution(2048);
  if (poly->Pixelize()) {
    std::cout << "Pixelization success!\nArea Check: Real: " <<
        poly->Area() << ", Pixelized: " << poly->PixelizedArea() << "\n";
  } else {
    std::cout << "Pixelization failed...\n";
  }
}

int main(int argc, char **argv) {

  std::string usage = "Usage: ";
  usage += argv[0];
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  // Check the FootprintBound class, specifically the derived class for making
  // circular footprints around a central point.
  if (FLAGS_all_footprint_tests || FLAGS_circle_footprint_tests)
    CircleFootprintTests();

  // Check the wedge FootprintBound class.
  if (FLAGS_all_footprint_tests || FLAGS_wedge_footprint_tests)
    WedgeFootprintTests();

  // Check the spherical polygon derived FootprintBound class.
  if (FLAGS_all_footprint_tests || FLAGS_polygon_footprint_tests)
    PolygonFootprintTests();

  return 0;
}
