#include <gtest/gtest.h>

#include "angular_bin-inl.h"

#include "pixel.h"
#include "point.h"

TEST(angular_bin, TestAngularBinConstructor) {
  // Test the default constructor, which should merely duplicate the input
  // angular bins, with no regionation or level set.
  double theta_min_deg = 0.01;
  double theta_max_deg = 1.0;
  double theta_mid_deg = 0.5 * (theta_min_deg + theta_max_deg);
  s2omp::angular_bin bin(theta_min_deg, theta_max_deg);

  // Check our bin bounds.
  ASSERT_DOUBLE_EQ(bin.theta_min(), theta_min_deg);
  ASSERT_DOUBLE_EQ(bin.theta_max(), theta_max_deg);
  ASSERT_TRUE(bin.is_within_bounds(theta_mid_deg));
  ASSERT_FALSE(bin.is_within_bounds(10.0 * theta_mid_deg));
  ASSERT_FALSE(bin.is_within_bounds(theta_min_deg + theta_max_deg));

  // Check the default values for the various counters and parameters.
  ASSERT_EQ(bin.n_region(), 0);
  ASSERT_EQ(bin.level(), -1);
  ASSERT_EQ(bin.pair_counts(), 0);
  ASSERT_DOUBLE_EQ(bin.pair_weight(), 0.0);
  ASSERT_DOUBLE_EQ(bin.pair_weight(s2omp::angular_bin::GAL_GAL), 0.0);
  ASSERT_DOUBLE_EQ(bin.pair_weight(s2omp::angular_bin::GAL_RAND), 0.0);
  ASSERT_DOUBLE_EQ(bin.pair_weight(s2omp::angular_bin::RAND_GAL), 0.0);
  ASSERT_DOUBLE_EQ(bin.pair_weight(s2omp::angular_bin::RAND_RAND), 0.0);
  ASSERT_DOUBLE_EQ(bin.pixel_wtheta(), 0.0);
  ASSERT_DOUBLE_EQ(bin.pixel_weight(), 0.0);

  // Now specify the number of regions.
  int n_regions = 10;
  bin = s2omp::angular_bin(theta_min_deg, theta_max_deg, n_regions);
  ASSERT_EQ(bin.n_region(), n_regions);
  ASSERT_TRUE(bin.regions_initialized());
  ASSERT_EQ(bin.level(), -1);
  // Default counter values should be 0 for all regions as well as region -1,
  // which should call the aggregate values over all regions.
  for (int k = -1; k < n_regions; k++) {
    ASSERT_EQ(bin.pair_counts(k), 0);
    ASSERT_DOUBLE_EQ(bin.pair_weight(k), 0.0);
    ASSERT_DOUBLE_EQ(bin.pair_weight(s2omp::angular_bin::GAL_GAL, k), 0.0);
    ASSERT_DOUBLE_EQ(bin.pair_weight(s2omp::angular_bin::GAL_RAND, k), 0.0);
    ASSERT_DOUBLE_EQ(bin.pair_weight(s2omp::angular_bin::RAND_GAL, k), 0.0);
    ASSERT_DOUBLE_EQ(bin.pair_weight(s2omp::angular_bin::RAND_RAND, k), 0.0);
    ASSERT_DOUBLE_EQ(bin.pixel_wtheta(k), 0.0);
    ASSERT_DOUBLE_EQ(bin.pixel_weight(k), 0.0);
  }
}

TEST(angular_bin, TestAngularBinLevelFinder) {
  // Test the routines for setting and finding the appropriate pixelization
  // level for our angular bin.
  double theta_min_deg = 0.001;
  double theta_max_deg = 10.0;
  s2omp::angular_bin bin(theta_min_deg, theta_max_deg);

  int level = 20;
  bin.set_level(level);
  ASSERT_EQ(bin.level(), level);

  // We don't want to hard code the levels used for each bin, but we can
  // verify that, if we have a set of bins that increase in scale monotonically
  // (logrithmically, in this case), then we expect to see a monotonic
  // decrease in the level for that bin.
  double bins_per_decade = 5.0;
  double unit_double = floor(log10(theta_min_deg)) * bins_per_decade;
  double theta_min = theta_min_deg;
  int current_level = s2omp::MAX_LEVEL;
  while (theta_min < theta_max_deg) {
    theta_min = pow(10.0, unit_double / bins_per_decade);
    double theta_max = pow(10.0, (unit_double + 1.0) / bins_per_decade);
    bin = s2omp::angular_bin(theta_min, theta_max);
    bin.find_level();
    std::cout << theta_min << " - " << theta_max << ": " << bin.level() << "\n";

    ASSERT_LE(bin.level(), current_level);

    current_level = bin.level();
    unit_double += 1.0;
  }
}
