#include <gtest/gtest.h>

//#include "point_test.cc"

//#include "pixel_test.cc"
//#include "field_pixel_test.cc"
//#include "tree_pixel_test.cc"

//#include "angular_bin_test.cc"
//#include "circle_bound_test.cc"
//#include "annulus_bound_test.cc"
#include "polygon_bound_test.cc"

//#include "pixel_union_test.cc"
//#include "field_union_test.cc"
//#include "tree_union_test.cc"

//#include "coverer_test.cc"
//#include "region_map_test.cc"
//#include "io_test.cc"

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
