#include <gtest/gtest.h>

#include "region_map.h"

#include "circle_bound.h"
#include "coverer.h"

bool VERBOSE = true;

TEST(region_map, TestRegionMapBasics) {
  // The default constructor without any arguments should produce an empty
  // region_map.
  s2omp::region_map mapper;

  ASSERT_TRUE(mapper.is_empty());
  ASSERT_EQ(mapper.level(), -1);
  ASSERT_EQ(mapper.n_region(), 0);

  // Even though the region_map is empty, we should still be able to call the
  // various methods for finding a region index without causing an error.  In
  // all cases, the return values should be INVALID_REGION_VALUE;
  s2omp::point p = s2omp::point::from_radec_deg(0.0, 0.0);

  ASSERT_EQ(mapper.find_region(p.to_pixel()), s2omp::INVALID_REGION_VALUE);
  ASSERT_EQ(mapper.find_region(p.to_pixel(0)), s2omp::INVALID_REGION_VALUE);
  ASSERT_EQ(mapper.find_region(p), s2omp::INVALID_REGION_VALUE);
  ASSERT_DOUBLE_EQ(mapper.get_area(0), -1.0);
}

double get_mean_area(s2omp::region_map& mapper) {
  double mean_area = 0.0;
  for (int k = 0; k < mapper.n_region(); k++) {
    mean_area += mapper.get_area(k) / mapper.n_region();
  }
  return mean_area;
}

double get_std_area(s2omp::region_map& mapper, double mean_area) {
  double std_area = 0.0;
  for (int k = 0; k < mapper.n_region(); k++) {
    std_area +=
        (mean_area - mapper.get_area(k)) * (mean_area - mapper.get_area(k));
  }
  return sqrt(std_area) / mapper.n_region();
}

void profile_region_map(std::string tag, s2omp::region_map& mapper) {
  std::cout << tag << " profile\n";
  std::cout << "\tn_region = " << mapper.n_region() << ", level = " <<
      mapper.level() << "\n";
  double mean_area = get_mean_area(mapper);
  double std_area = get_std_area(mapper, mean_area);
  std::cout << "\tMean area = " << mean_area << " +- " <<
      sqrt(std_area) / mapper.n_region() << " sq. deg.\n";
}

bool valid_region(int region) {
  return region != s2omp::INVALID_REGION_VALUE;
}

TEST(region_map, TestRegionMapRegionation) {
  // Test the basics of sub-dividing a bound_interface object into roughly
  // equal area pieces.

  // Start with a circle_bound
  // Start by constructing a simple circle_bound.
  s2omp::point axis = s2omp::point::from_radec_deg(0.0, 0.0);
  double theta_deg = 10.0;
  s2omp::circle_bound* bound =
      s2omp::circle_bound::from_radius(axis, theta_deg);

  // Construct a region_map for this bound with the default parameters
  s2omp::region_map mapper;
  int n_region = 10;
  uint n_region_actual = mapper.init(*bound, n_region);
  ASSERT_GE(mapper.level(), 0);
  ASSERT_LE(mapper.level(), s2omp::MAX_LEVEL);
  ASSERT_EQ(mapper.n_region(), n_region_actual);
  EXPECT_EQ(n_region_actual, n_region);
  if (VERBOSE) profile_region_map("BasicRegionMap", mapper);

  // Construct a simple covering of the bound at the same level and verify
  // that all covering pixels are within the region_map.
  s2omp::pixel_vector covering;
  s2omp::coverer::get_simple_covering(*bound, mapper.level(), &covering);
  for (s2omp::pixel_iterator iter = covering.begin();
      iter != covering.end(); ++iter) {
    ASSERT_TRUE(valid_region(mapper.find_region(*iter)));
    ASSERT_TRUE(valid_region(mapper.find_region(iter->center_point())));
    // Children of covering pixels should produce a valid region index, while
    // parents should not.
    ASSERT_TRUE(valid_region(mapper.find_region(iter->child_begin())));
    ASSERT_FALSE(valid_region(mapper.find_region(iter->parent())));
  }

  // In addition to the center points of the covering pixels, we should be able
  // to find a region index for any point within the starting bound.
  int n_random = 1000;
  for (int k = 0; k < n_random; k++) {
    ASSERT_TRUE(valid_region(mapper.find_region(bound->create_random_point())));
  }

  // Now produce a region_map at a coarser resolution.  This should still be a
  // valid regionation, but the chances that we have regions with wildly
  // different areas are greater.
  int level = mapper.level() - 2;
  s2omp::region_map coarse_mapper;
  n_region_actual = coarse_mapper.init(*bound, n_region, level);
  ASSERT_EQ(coarse_mapper.level(), level);
  ASSERT_EQ(coarse_mapper.n_region(), n_region_actual);
  EXPECT_EQ(n_region_actual, n_region);
  if (VERBOSE) profile_region_map("CoarseRegionMap", coarse_mapper);

  // Verify with a simple covering.
  s2omp::coverer::get_simple_covering(*bound, coarse_mapper.level(), &covering);
  for (s2omp::pixel_iterator iter = covering.begin();
      iter != covering.end(); ++iter) {
    ASSERT_TRUE(valid_region(coarse_mapper.find_region(*iter)));
    ASSERT_TRUE(valid_region(coarse_mapper.find_region(iter->center_point())));
    ASSERT_TRUE(valid_region(coarse_mapper.find_region(iter->child_begin())));
    ASSERT_FALSE(valid_region(coarse_mapper.find_region(iter->parent())));
  }

  // And for good measure, generate a region map at finer resolution.
  level = mapper.level() + 2;
  s2omp::region_map fine_mapper;
  n_region_actual = fine_mapper.init(*bound, n_region, level);
  ASSERT_EQ(fine_mapper.level(), level);
  ASSERT_EQ(fine_mapper.n_region(), n_region_actual);
  EXPECT_EQ(n_region_actual, n_region);
  if (VERBOSE) profile_region_map("FineRegionMap", fine_mapper);
  // Verify that, by increasing the resolution of our region_map, we've reduced
  // the variance in area between regions.
  ASSERT_LT(get_std_area(fine_mapper, get_mean_area(fine_mapper)),
      get_std_area(mapper, get_mean_area(mapper)));
}
