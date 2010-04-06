## STOMP source directory

# Turn off README, NEWS, and other GNU files.
AUTOMAKE_OPTIONS = foreign

# Common flags.
CXXFLAGS = @CXXFLAGS@ -Wall
LDFLAGS = @LDFLAGS@ @GFLAGS_LIB@ -lgflags -lstomp
INCLUDES = -I@top_srcdir@/stomp/ @GFLAGS_INCLUDE@

h_sources = MersenneTwister.h stomp_angular_bin.h stomp_angular_coordinate.h stomp_angular_correlation.h stomp_base_map.h stomp_core.h stomp_geometry.h stomp_map.h stomp_pixel.h stomp_scalar_map.h stomp_scalar_pixel.h stomp_tree_map.h stomp_tree_pixel.h stomp_util.h
cc_sources = stomp_angular_bin.cc stomp_angular_coordinate.cc stomp_angular_correlation.cc stomp_base_map.cc stomp_core.cc stomp_geometry.cc stomp_map.cc stomp_pixel.cc stomp_scalar_map.cc stomp_scalar_pixel.cc stomp_tree_map.cc stomp_tree_pixel.cc stomp_util.cc

library_includedir=$(includedir)/$(GENERIC_LIBRARY_NAME)-$(GENERIC_API_VERSION)/$(GENERIC_LIBRARY_NAME)
library_include_HEADERS = $(h_sources)

toplevel_includedir=$(includedir)
toplevel_include_HEADERS = stomp.h

INCLUDES = -I$(top_srcdir)

lib_LTLIBRARIES= libstomp-1.0.la
libstomp_1_0_la_SOURCES= $(h_sources) $(cc_sources)
libstomp_1_0_la_LDFLAGS= -version-info $(GENERIC_LIBRARY_VERSION) -release
$(GENERIC_RELEASE)

# Test programs
check_PROGRAMS = stomp_angular_coordinate_test stomp_angular_correlation_test stomp_core_test stomp_geometry_test stomp_map_test stomp_pixel_test stomp_scalar_map_test stomp_scalar_pixel_test stomp_tree_map_test stomp_tree_pixel_test stomp_unit_test stomp_util_test

bin_PROGRAMS = stomp_unit_test
stomp_unit_test_SOURCES = stomp_angular_coordinate_test.cc stomp_angular_correlation_test.cc stomp_core_test.cc stomp_geometry_test.cc stomp_map_test.cc stomp_pixel_test.cc stomp_scalar_map_test.cc stomp_scalar_pixel_test.cc stomp_tree_map_test.cc stomp_tree_pixel_test.cc stomp_util_test.cc stomp_unit_test.cc
stomp_unit_test_LDADD = libstomp-1.0.la

# Test programs run automatically by 'make check'
TESTS = $(check_PROGRAMS)

