#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>

// Define our command-line flags.
DEFINE_bool(all_tests, false, "Run all unit tests.");

int main(int argc, char **argv) {
  void CoreUnitTests(bool run_all_tests);
  void AngularCoordinateUnitTests(bool run_all_tests);
  void AngularCorrelationUnitTests(bool run_all_tests);
  void PixelUnitTests(bool run_all_tests);
  void ScalarPixelUnitTests(bool run_all_tests);
  void TreePixelUnitTests(bool run_all_tests);
  void IndexedTreePixelUnitTests(bool run_all_tests);
  void MapUnitTests(bool run_all_tests);
  void ScalarMapUnitTests(bool run_all_tests);
  void TreeMapUnitTests(bool run_all_tests);
  void IndexedTreeMapUnitTests(bool run_all_tests);
  void GeometryUnitTests(bool run_all_tests);
  void UtilUnitTests(bool run_all_tests);

  std::string usage = "Usage: ";
  usage += argv[0];
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  // Now we call each of the unit test suites, which will run the requested
  // unit tests.

  // The Core class
  CoreUnitTests(FLAGS_all_tests);

  // The AngularCoordinate class and derived classes
  AngularCoordinateUnitTests(FLAGS_all_tests);

  // The AngularCorrelation and AngularBin classes
  AngularCorrelationUnitTests(FLAGS_all_tests);

  // The Pixel class
  PixelUnitTests(FLAGS_all_tests);

  // The ScalarPixel class
  ScalarPixelUnitTests(FLAGS_all_tests);

  // The TreePixel class
  TreePixelUnitTests(FLAGS_all_tests);

  // The IndexedTreePixel class
  IndexedTreePixelUnitTests(FLAGS_all_tests);

  // The Map class
  MapUnitTests(FLAGS_all_tests);

  // The ScalarMap class
  ScalarMapUnitTests(FLAGS_all_tests);

  // The TreeMap class
  TreeMapUnitTests(FLAGS_all_tests);

  // The IndexedTreeMap class
  IndexedTreeMapUnitTests(FLAGS_all_tests);

  // The GeometricBound class and derivatives
  GeometryUnitTests(FLAGS_all_tests);

  // The utility classes
  UtilUnitTests(FLAGS_all_tests);

  return 0;
}


