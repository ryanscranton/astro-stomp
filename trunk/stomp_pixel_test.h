// Define our command line flags here so we can use these flags later.
DEFINE_bool(all_pixel_tests, false, "Run all class unit tests.");
DEFINE_bool(pixel_basic_tests, false, "Run Pixel resolution tests");
DEFINE_bool(pixel_stripe_tests, false, "Run Pixel stripe tests");
DEFINE_bool(pixel_xy_tests, false, "Run Pixel XY tests");
DEFINE_bool(pixel_bound_tests, false, "Run Pixel bound tests");
DEFINE_bool(pixel_within_radius_tests, false, "Run Pixel WithinRadius tests");
DEFINE_bool(pixel_annulus_intersection_tests, false,
            "Run Pixel AnnulusIntersection tests");
