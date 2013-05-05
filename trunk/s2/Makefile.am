## Forked S2 source directory

# Turn off README, NEWS, and other GNU files.
AUTOMAKE_OPTIONS = foreign

# Common flags.
CXXFLAGS = @CXXFLAGS@ -Wall -fPIC -O -I. -DARCH_K8 -Wno-deprecated -DHASH_NAMESPACE=__gnu_cxx
LDFLAGS = @LDFLAGS@ -arch x86_64 -bundle -flat_namespace -undefined suppress @GFLAGS_LIB@ @ZLIB_LIB@
INCLUDES = -I@top_srcdir@/s2/ @GFLAGS_INCLUDE@ @ZLIB_INCLUDE@

h_sources = base/int128.h base/logging.h base/stringprintf.h util/coding/coder.h util/coding/varint.h util/math/mathutil.h util/math/mathlimits.h s1angle.h s2.h s2cellid.h s2latlng.h s1interval.h s2cap.h s2cell.h s2cellunion.h s2edgeindex.h s2edgeutil.h s2latlngrect.h s2loop.h s2pointregion.h s2polygon.h s2polygonbuilder.h s2polyline.h s2r2rect.h s2region.h s2regioncoverer.h s2regionintersection.h s2regionunion.h

cc_sources = base/int128.cc base/logging.cc base/stringprintf.cc util/coding/coder.cc util/coding/varint.cc util/math/mathutil.cc util/math/mathlimits.cc s1angle.cc s2.cc s2cellid.cc s2latlng.cc s1interval.cc s2cap.cc s2cell.cc s2cellunion.cc s2edgeindex.cc s2edgeutil.cc s2latlngrect.cc s2loop.cc s2pointregion.cc s2polygon.cc s2polygonbuilder.cc s2polyline.cc s2r2rect.cc s2region.cc s2regioncoverer.cc s2regionintersection.cc s2regionunion.cc

library_includedir=$(includedir)/s2
library_include_HEADERS = s1angle.h s2.h s2cellid.h s2latlng.h s1interval.h s2cap.h s2cell.h s2cellunion.h s2edgeindex.h s2edgeutil.h s2latlngrect.h s2loop.h s2pointregion.h s2polygon.h s2polygonbuilder.h s2polyline.h s2r2rect.h s2region.h s2regioncoverer.h s2regionintersection.h s2regionunion.h

base_includedir=$(includedir)/base
base_include_HEADERS = base/int128.h base/logging.h base/stringprintf.h

coding_includedir=$(includedir)/util/coding
coding_include_HEADERS = util/coding/coder.h util/coding/varint.h

math_includedir=$(includedir)/util/coding
math_include_HEADERS = util/math/mathutil.h util/math/mathlimits.h

lib_LTLIBRARIES= libs2.la
libs2_la_SOURCES= $(h_sources) $(cc_sources)
libs2_la_LDFLAGS= -version-info $(GENERIC_LIBRARY_VERSION) -release $(GENERIC_RELEASE)

