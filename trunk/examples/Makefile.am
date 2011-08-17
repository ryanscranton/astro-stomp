## STOMP source directory

# Turn off README, NEWS, and other GNU files.
AUTOMAKE_OPTIONS = foreign

# Common flags.
CXXFLAGS = @CXXFLAGS@ -Wall
LDFLAGS = @LDFLAGS@ -lstomp @GFLAGS_LIB@
INCLUDES = -I@top_srcdir@/stomp/ @GFLAGS_INCLUDE@
LIBS = -L@top_srcdir@/stomp/ @GFLAGS_LIB@

exampledir = @top_srcdir@/examples

example_PROGRAMS = stomp_adapt_map stomp_circle_coverage stomp_combine_maps \
                   stomp_galaxy_autocorrelation \
                   stomp_galaxy_autocorrelation_jack stomp_local_density \
                   stomp_galaxy_radial_autocorrelation_jack \
                   stomp_pixelize_circles stomp_stripe_filter \
                   stomp_sysmap_density stomp_sysmap_filter \
                   stomp_sysmap_generator stomp_sys_xcorr

stomp_adapt_map_SOURCES = stomp_adapt_map.cc

stomp_circle_coverage_SOURCES = stomp_circle_coverage.cc

stomp_combine_maps_SOURCES = stomp_combine_maps.cc

stomp_galaxy_autocorrelation_SOURCES = stomp_galaxy_autocorrelation.cc

stomp_galaxy_autocorrelation_jack_SOURCES = stomp_galaxy_autocorrelation_jack.cc

stomp_galaxy_radial_autocorrelation_jack_SOURCES = stomp_galaxy_radial_autocorrelation_jack.cc

stomp_local_density_SOURCES = stomp_local_density.cc

stomp_pixelize_circles_SOURCES = stomp_pixelize_circles.cc

stomp_projected_density_SOURCES = stomp_projected_density.cc

stomp_stripe_filter_SOURCES = stomp_stripe_filter.cc

stomp_sysmap_density_SOURCES = stomp_sysmap_density.cc

stomp_sysmap_filter_SOURCES = stomp_sysmap_filter.cc

stomp_sysmap_generator_SOURCES = stomp_sysmap_generator.cc

stomp_sys_xcorr_SOURCES = stomp_sys_xcorr.cc
