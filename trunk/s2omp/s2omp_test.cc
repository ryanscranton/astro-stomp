#include <gtest/gtest.h>

#include "point_test.cc"

#include "pixel_test.cc"

#include "angular_bin_test.cc"
#include "circle_bound_test.cc"
#include "annulus_bound_test.cc"

#include "pixel_union_test.cc"

#include "coverer_test.cc"
#include "region_map_test.cc"

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
