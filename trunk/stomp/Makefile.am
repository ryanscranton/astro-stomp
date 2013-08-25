## STOMP source directory

# Turn off README, NEWS, and other GNU files.
AUTOMAKE_OPTIONS = foreign

# Common flags.
CXXFLAGS = @CXXFLAGS@ -Wall -std=c++0x -stdlib=libc++
LDFLAGS = @LDFLAGS@ @GFLAGS_LIB@
INCLUDES = -I@top_srcdir@/stomp/ @GFLAGS_INCLUDE@

h_sources = MersenneTwister.h stomp_angular_bin.h stomp_angular_coordinate.h stomp_angular_correlation.h stomp_base_map.h stomp_core.h stomp_geometry.h stomp_map.h stomp_pixel.h stomp_scalar_map.h stomp_scalar_pixel.h stomp_tree_map.h stomp_tree_pixel.h stomp_util.h stomp_itree_pixel.h stomp_itree_map.h stomp_radial_bin.h stomp_radial_correlation.h
cc_sources = stomp_angular_bin.cc stomp_angular_coordinate.cc stomp_angular_correlation.cc stomp_base_map.cc stomp_core.cc stomp_geometry.cc stomp_map.cc stomp_pixel.cc stomp_scalar_map.cc stomp_scalar_pixel.cc stomp_tree_map.cc stomp_tree_pixel.cc stomp_util.cc stomp_itree_pixel.cc stomp_itree_map.cc stomp_radial_bin.cc stomp_radial_correlation.cc

library_includedir=$(includedir)/$(GENERIC_LIBRARY_NAME)/
library_include_HEADERS = $(h_sources)

toplevel_includedir=$(includedir)
toplevel_include_HEADERS = stomp.h

lib_LTLIBRARIES= libstomp.la
libstomp_la_SOURCES= $(h_sources) $(cc_sources)
libstomp_la_LDFLAGS= -version-info $(GENERIC_LIBRARY_VERSION) -release $(GENERIC_RELEASE)

check_PROGRAMS = stomp_unit_test
stomp_unit_test_SOURCES = stomp_angular_coordinate_test.cc stomp_angular_correlation_test.cc stomp_core_test.cc stomp_geometry_test.cc stomp_map_test.cc stomp_pixel_test.cc stomp_scalar_map_test.cc stomp_scalar_pixel_test.cc stomp_tree_map_test.cc stomp_itree_map_test.cc stomp_tree_pixel_test.cc stomp_itree_pixel_test.cc stomp_util_test.cc stomp_unit_test.cc
stomp_unit_test_LDADD = libstomp.la

# Test programs run automatically by 'make check'
TESTS = $(check_PROGRAMS)

