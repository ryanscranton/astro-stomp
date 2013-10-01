## STOMP source directory

# Turn off README, NEWS, and other GNU files.
AUTOMAKE_OPTIONS = foreign

# Common flags.
CXXFLAGS = @CXXFLAGS@ -Wall -std=c++0x -stdlib=libc++ -Wno-deprecated -DHASH_NAMESPACE=__gnu_cxx -Wgnu
LDFLAGS = @LDFLAGS@ @GFLAGS_LIB@ -L@top_srcdir@/s2 -ls2 -lprotobuf
INCLUDES = -I@top_srcdir@/s2omp -I@top_srcdir@/stomp -I@top_srcdir@/s2 @GFLAGS_INCLUDE@

h_sources = point.h pixel.h bound_interface.h coverer.h circle_bound.h angular_bin-inl.h annulus_bound.h region_map.h util.h cosmo_point-inl.h indexed_point-inl.h pixel_union.h tree_pixel.h tree_union.h field_pixel-inl.h field_union.h latlon_bound.h io.h s2omp.pb.h
cc_sources = point.cc pixel.cc bound_interface.cc coverer.cc circle_bound.cc annulus_bound.cc region_map.cc pixel_union.cc tree_pixel.cc tree_union.cc field_union.cc latlon_bound.cc io.cc s2omp.pb.cc

library_includedir=$(includedir)/s2omp
library_include_HEADERS = $(h_sources)

toplevel_includedir=$(includedir)
toplevel_include_HEADERS = s2omp.h

lib_LTLIBRARIES= libs2omp.la
libs2omp_la_SOURCES= $(h_sources) $(cc_sources)
libs2omp_la_LDFLAGS= -version-info $(GENERIC_LIBRARY_VERSION) -release $(GENERIC_RELEASE)

check_PROGRAMS = s2omp_test
s2omp_test_SOURCES = s2omp_test.cc
s2omp_test_LDADD = libs2omp.la -lgtest -lgtest_main -lpthread

# Test programs run automatically by 'make check'
TESTS = $(check_PROGRAMS)

