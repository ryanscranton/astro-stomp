#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <stomp.h>
#include <gflags/gflags.h>

DEFINE_string(map_file, "",
              "Name of the ASCII file containing the StompMap geometry");
DEFINE_string(galaxy_file, "",
              "Name of the ASCII file containing Galaxies");
DEFINE_bool(galaxy_radec, false, "Galaxy coordinates are in RA-DEC");
DEFINE_string(output_file, "Wradial_test",
              "Output file");
DEFINE_double(r_min, 1, "Minimum radial scale (in Mpc)");
DEFINE_double(r_max, 10.0, "Maximum radial scale (in Mpc)");
DEFINE_double(mag_min, 16.0, "Minimum acceptable galaxy magnitude");
DEFINE_double(mag_max, 17.0, "Maximum acceptable galaxy magnitude");
DEFINE_double(z_min, 0.0, "Minimum acceptable galaxy redshift");
DEFINE_double(z_max, 5.01, "Maximum acceptable galaxy redshift");
DEFINE_double(prob_min, 0.2, "Minimum acceptable galaxy likelihood");
DEFINE_double(prob_max, 1.00001, "Maximum acceptable galaxy likelihood");
DEFINE_int32(n_bins_per_decade, 5, "Number of angular bins per decade.");
DEFINE_int32(n_random, 1,
	     "Integer number of random points per galaxy to use.");
DEFINE_int32(n_jackknife, 0,
	     "Number of jack-knife samples to use. Defaults to 2*angular bins");
DEFINE_bool(single_index, false, "Use older single-index file format.");
DEFINE_bool(no_weight, false, "Map file is missing weight column.");
DEFINE_int32(maximum_resolution, -1,
	     "Maximum resolution to use for pixel-based estimator");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --map_file=<StompMap ASCII>";
  usage += " --galaxy_file=<Galaxy catalog ASCII>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_map_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_galaxy_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  // First, we read our STOMP map into a map object.  There are a couple
  // permutations based on the various map formats that are out there: with
  // or without a weight column or in the single index or double index format.
  Stomp::Map* stomp_map = new Stomp::Map(FLAGS_map_file, !FLAGS_single_index,
					 !FLAGS_no_weight);
  std::cout << "Read map from " << FLAGS_map_file << "; total area: " <<
    stomp_map->Area() << " sq. deg.\n";

  // Now we read in our galaxy data file.  The expected format is
  //  THETA  PHI  WEIGHT  MAGNITUDE REDSHIFT
  // where the WEIGHT column is the likelihood that the object is a galaxy
  // and MAGNITUDE is the apparent magnitude in a given filter.  We filter all
  // of the objects against the map, tossing out any objects that aren't in the
  // map.
  Stomp::CosmoVector galaxy;
  std::ifstream galaxy_file(FLAGS_galaxy_file.c_str());
  double theta, phi, prob, mag, redshift;
  uint32_t n_galaxy = 0;

  Stomp::AngularCoordinate::Sphere galaxy_sphere =
    Stomp::AngularCoordinate::Survey;
  if (FLAGS_galaxy_radec) galaxy_sphere = Stomp::AngularCoordinate::Equatorial;
  while (!galaxy_file.eof()) {
    galaxy_file >> theta >> phi >> prob >> mag >> redshift;
    Stomp::CosmoCoordinate tmp_ang(theta, phi, prob, redshift,
				   galaxy_sphere);

    if ((prob >= FLAGS_prob_min) && (prob <= FLAGS_prob_max) && 
	(mag >= FLAGS_mag_min) && (mag <= FLAGS_mag_max) && 
	(redshift >= FLAGS_z_min) && (redshift <= FLAGS_z_max) &&
	stomp_map->Contains(tmp_ang) && (tmp_ang.Weight() > 0.2))
      galaxy.push_back(tmp_ang);
    n_galaxy++;
  }
  galaxy_file.close();

  std::cout << "Read " << n_galaxy << " galaxies; kept " <<
    galaxy.size() << "\n";
  n_galaxy = galaxy.size();
  galaxy.resize(n_galaxy);

  // Now, we set up the object that will contain the measurement results.  The
  // correlation object is a essentially a container for radial bin objects
  // which have a given radial range (all object or pixel pairs separated by
  // 1 < Mpc/h < 10 degrees, for instance).  In addition, the constructor
  // for these objects will work out, based on the radial bin size, which
  // StompMap resolution would be appropriate for calculating the radial
  // correlation on that scale.
  Stomp::RadialCorrelation wradial(FLAGS_r_min, FLAGS_r_max,
				  FLAGS_n_bins_per_decade);
  //Stomp::RadialCorrelation wtheta(FLAGS_theta_min, FLAGS_theta_max,
  //				  FLAGS_n_bins_per_decade);

  wradial.FindAutoCorrelationWithRegions(*stomp_map, galaxy,
					 static_cast<uint8_t>(FLAGS_n_random),
					 static_cast<uint16_t>(FLAGS_n_jackknife));

  wradial.Write(FLAGS_output_file);
  
  return 0;
}
