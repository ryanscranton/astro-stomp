#include "stomp_core_test.h"
#include "stomp_angular_coordinate_test.h"
#include "stomp_angular_correlation_test.h"
#include "stomp_pixel_test.h"
#include "stomp_scalar_pixel_test.h"
#include "stomp_tree_pixel_test.h"
#include "stomp_map_test.h"
#include "stomp_scalar_map_test.h"
#include "stomp_tree_map_test.h"
#include "stomp_footprint_test.h"
#include "stomp_util_test.h"
#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>

// Define our command-line flags.
DEFINE_bool(all_tests, false, "Run all unit tests.");

int main(int argc, char **argv) {
  void ConstantsTests();
  void MSBTests();
  void CosmologyTests();
  void AngularCoordinateTests();
  void AngularCoordinatePositionAngleTests();
  void AngularCoordinateRotationTests();
  void WeightedAngularCoordinateTests();
  void PixelResolutionTests();
  void PixelStripeTests();
  void PixelXYTests();
  void PixelBoundTests();
  void PixelWithinRadiusTests();
  void PixelAnnulusIntersectionTests();
  void ScalarPixelBasicTests();
  void TreePixelBasicTests();
  void TreePixelPairTests();
  void TreePixelCoverageTests();
  void TreePixelFieldPairTests();
  void TreePixelNeighborTests();
  void MapBasicTests();
  void MapWriteTests();
  void MapReadTests();
  void MapCoverTests();
  void MapIteratorTests();
  void MapLocationTests();
  void MapUnmaskedFractionTests();
  void MapContainsTests();
  void MapRandomPointsTests();
  void MapMultiMapTests();
  void MapRegionTests();
  void MapSoftenTests();
  void ScalarMapBasicTests();
  void ScalarMapLocalTests();
  void ScalarMapResamplingTests();
  void ScalarMapRegionTests();
  void ScalarMapAutoCorrelationTests();
  void TreeMapBasicTests();
  void TreeMapPairTests();
  void TreeMapAreaTests();
  void TreeMapRegionTests();
  void TreeMapFieldPairTests();
  void TreeMapNeighborTests();
  void AngularBinningTests();
  void CircleFootprintTests();
  void WedgeFootprintTests();
  void PolygonFootprintTests();

  std::string usage = "Usage: ";
  usage += argv[0];
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  // First things first.  Let's check to make sure that our constants are
  // all in place.
  if (FLAGS_all_tests || FLAGS_constants_tests) ConstantsTests();

  // Now, we check our MostSignificantBit method to make sure it's working.
  if (FLAGS_all_tests || FLAGS_msb_tests) MSBTests();

  // Now, we check our static Cosmology class to make sure it's functioning.
  if (FLAGS_all_tests || FLAGS_cosmology_tests) CosmologyTests();

  // Check that AngularCoordinate transforms from one coordinate system to
  // another work properly and that we can recover the same pixel indices from
  // transformed coordinates.
  if (FLAGS_all_tests || FLAGS_all_angular_coordinate_tests ||
      FLAGS_angular_coordinate_tests) AngularCoordinateTests();

  // Check that the AngularCoordinate position angles are being calculated
  // correctly.
  if (FLAGS_all_tests || FLAGS_all_angular_coordinate_tests ||
      FLAGS_angular_coordinate_position_angle_tests)
    AngularCoordinatePositionAngleTests();

  // Check that the AngularCoordinate rotations are being calculated
  // correctly.
  if (FLAGS_all_tests || FLAGS_all_angular_coordinate_tests ||
      FLAGS_angular_coordinate_rotation_tests)
    AngularCoordinateRotationTests();

  // Check WeightedAngularCoordinate extensions to the basic AngularCoordinate.
  if (FLAGS_all_tests || FLAGS_all_angular_coordinate_tests ||
      FLAGS_weighted_angular_coordinate_tests)
    WeightedAngularCoordinateTests();

  // Check that the hierarchical resolution scaling routines work properly.
  if (FLAGS_all_tests || FLAGS_all_pixel_tests ||
      FLAGS_pixel_resolution_tests) PixelResolutionTests();

  // Check that the routines for generating the X-Y bounds of a given
  // region work correctly.
  if (FLAGS_all_tests || FLAGS_all_pixel_tests ||
      FLAGS_pixel_xy_tests) PixelXYTests();

  // Quick checks for the various routines to return the pixel bounds in the
  // various coordinate systems.
  if (FLAGS_all_tests || FLAGS_all_pixel_tests ||
      FLAGS_pixel_bound_tests) PixelBoundTests();

  // Check the routines for taking the X-Y bounds and returning a list of pixels
  // with the specified radius.
  if (FLAGS_all_tests || FLAGS_all_pixel_tests ||
      FLAGS_within_radius_tests) WithinRadiusTests();

  // Check the routines for determining whether or not annuli intersect pixels.
  if (FLAGS_all_tests || FLAGS_all_pixel_tests ||
      FLAGS_annulus_intersection_tests) AnnulusIntersectionTests();

  // Check the basic inheritance of the Stomp::ScalarPixel class from
  // Stomp::Pixel.
  if (FLAGS_all_tests || FLAGS_all_scalar_pixel_tests ||
      FLAGS_scalar_pixel_tests) ScalarPixelTests();

  // Check that the Stomp::TreePixel class is able to add points and
  // automatically generate sub-pixels.
  if (FLAGS_all_tests || FLAGS_all_tree_pixel_tests ||
      FLAGS_tree_pixel_tests) TreePixelTests();

  // Check the Stomp::TreePixel pair-finding routines.
  if (FLAGS_all_tests || FLAGS_all_tree_pixel_tests ||
      FLAGS_tree_pixel_pair_tests) TreePixelPairTests();

  // Check the Stomp::TreePixel Coverage method works as advertised.
  if (FLAGS_all_tests || FLAGS_all_tree_pixel_tests ||
      FLAGS_tree_pixel_coverage_tests) TreePixelCoverageTests();

  // Checking pair finding routines with Field values.
  if (FLAGS_all_tests || FLAGS_all_tree_pixel_tests ||
      FLAGS_tree_pixel_field_pair_tests) TreePixelFieldPairTests();

  // Checking nearest neighbor finding routines.
  if (FLAGS_all_tests || FLAGS_all_tree_pixel_tests ||
      FLAGS_tree_pixel_neighbor_tests) TreePixelNeighborTests();

  // Check the basic routines for generating Stomp::Map instances from a list
  // of contiguous pixels.
  if (FLAGS_all_tests || FLAGS_all_map_tests ||
      FLAGS_map_initialization_tests) MapInitializationTests();

  // Check the routines for writing Stomp::Map's to file.
  if (FLAGS_all_tests || FLAGS_all_map_tests ||
      FLAGS_map_write_tests) MapWriteTests();

  // Check the routines for reading a Stomp::Map from a simple ASCII file.
  if (FLAGS_all_tests || FLAGS_all_map_tests ||
      FLAGS_map_read_tests) MapReadTests();

  // Check the routines for finding the Stomp::Map Coverage and Covering.
  if (FLAGS_all_tests || FLAGS_all_map_tests ||
      FLAGS_map_cover_tests) MapCoverTests();

  // Check the routines for iterating through a Stomp::Map's pixels using the
  // Begin, End and Iterate methods.
  if (FLAGS_all_tests || FLAGS_all_map_tests ||
      FLAGS_map_iterator_tests) MapIteratorTests();

  // Check the Stomp::Map methods for checking locations against the area
  // covered by a Stomp::Map.
  if (FLAGS_all_tests || FLAGS_all_map_tests ||
      FLAGS_map_location_tests) MapLocationTests();

  // Check the routines for finding the area of a pixel covered by a Stomp::Map.
  if (FLAGS_all_tests || FLAGS_all_map_tests ||
      FLAGS_map_unmasked_fraction_tests) MapUnmaskedFractionTests();

  // Check the routines for testing whether one Stomp::Map is contained in
  // another Stomp::Map.
  if (FLAGS_all_tests || FLAGS_all_map_tests ||
      FLAGS_map_contains_tests) MapContainsTests();

  // Check the routine for generating random points within a Stomp::Map's area.
  if (FLAGS_all_tests || FLAGS_all_map_tests ||
      FLAGS_map_random_points_tests) MapRandomPointsTests();

  // Check the different ways of combining Stomp::Map instances.
  if (FLAGS_all_tests || FLAGS_all_map_tests ||
      FLAGS_map_multimap_tests) MapMultiMapTests();

  // Check the routines for breaking a Stomp::Map into sub-regions.
  if (FLAGS_all_tests || FLAGS_all_map_tests ||
      FLAGS_map_region_tests) MapRegionTests();

  // Check the routines for softening the maximum resolution of the
  // Map and cutting the map based on the Weight.
  if (FLAGS_all_tests || FLAGS_all_map_tests ||
      FLAGS_map_soften_tests) MapSoftenTests();

  // Check the basic routines for generating a Stomp::ScalarMap from an input
  // Stomp::Map.
  if (FLAGS_all_tests || FLAGS_all_scalar_map_tests ||
      FLAGS_scalar_map_tests) ScalarMapTests();

  // Check the StomScalarMap methods for finding the area and density of the
  // map within a given pixel.
  if (FLAGS_all_tests || FLAGS_all_scalar_map_tests ||
      FLAGS_scalar_map_local_tests) ScalarMapLocalTests();

  // Check the Stomp::ScalarMap methods for creating new, coarser resolution
  // Stomp::ScalarMaps from an initial high-resolution version.
  if (FLAGS_all_tests || FLAGS_all_scalar_map_tests ||
      FLAGS_scalar_map_resampling_tests) ScalarMapResamplingTests();

  // Check the routines for splitting up the area of a Stomp::ScalarMap into
  // roughly equal-area regions.
  if (FLAGS_all_tests || FLAGS_all_scalar_map_tests ||
      FLAGS_scalar_map_region_tests) ScalarMapRegionTests();

  // Check the auto-correlation methods in the Stomp::ScalarMap class.
  if (FLAGS_all_tests || FLAGS_all_scalar_map_tests ||
      FLAGS_scalar_map_autocorrelation_tests)
    StompScalarMapAutoCorrelationTests();

  // Check that the Stomp::TreeMap class is able to add points and
  // automatically generate sub-pixels.
  if (FLAGS_all_tests || FLAGS_all_tree_map_tests ||
      FLAGS_tree_map_tests) TreeMapTests();

  // Check the Stomp::TreeMap pair-finding routines.
  if (FLAGS_all_tests || FLAGS_all_tree_map_tests ||
      FLAGS_tree_map_pair_tests) TreeMapPairTests();

  // Check that the Stomp::TreeMap class area calculations improve as more
  // points are added to the map.
  if (FLAGS_all_tests || FLAGS_all_tree_map_tests ||
      FLAGS_tree_map_area_tests) TreeMapAreaTests();

  // Check that the Stomp::TreeMap class is able to regionate properly.
  if (FLAGS_all_tests || FLAGS_all_tree_map_tests ||
      FLAGS_tree_map_region_tests) TreeMapRegionTests();

  // Checking pair finding routines with Field values.
  if (FLAGS_all_tests || FLAGS_all_tree_map_tests ||
      FLAGS_tree_map_field_pair_tests) TreeMapFieldPairTests();

  // Checking nearest neighbor routines.
  if (FLAGS_all_tests || FLAGS_all_tree_map_tests ||
      FLAGS_tree_map_neighbor_tests) TreeMapNeighborTests();

  // Check the routines related to the AngularBin and AngularCorrelation
  // classes.
  if (FLAGS_all_tests || FLAGS_angular_bin_tests) AngularBinningTests();

  // Check the FootprintBound class, specifically the derived class for making
  // circular footprints around a central point.
  if (FLAGS_all_tests || FLAGS_all_footprint_tests ||
      FLAGS_circle_footprint_tests) CircleFootprintTests();

  // Check the wedge FootprintBound class.
  if (FLAGS_all_tests || FLAGS_all_footprint_tests ||
      FLAGS_wedge_footprint_tests) WedgeFootprintTests();

  // Check the spherical polygon derived FootprintBound class.
  if (FLAGS_all_tests || FLAGS_all_footprint_tests ||
      FLAGS_polygon_footprint_tests) PolygonFootprintTests();

  return 0;
}


