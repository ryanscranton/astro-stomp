// Define our command line flags here so we can use these flags later.
DEFINE_bool(all_scalar_map_tests, false, "Run all class unit tests.");
DEFINE_bool(scalar_map_basic_tests, false, "Run ScalarMap basic tests");
DEFINE_bool(scalar_map_local_tests, false, "Run ScalarMap local tests");
DEFINE_bool(scalar_map_resampling_tests, false,
            "Run ScalarMap resampling tests");
DEFINE_bool(scalar_map_region_tests, false, "Run ScalarMap region tests");
DEFINE_bool(scalar_map_autocorrelation_tests, false,
            "Run ScalarMap auto-correlation tests");
