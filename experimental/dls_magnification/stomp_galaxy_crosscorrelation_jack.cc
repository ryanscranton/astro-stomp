#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <stomp.h>
#include <gflags/gflags.h>

// Define our command-line flags.
DEFINE_string(map_file, "",
              "Name of the ASCII file containing the StompMap geometry");
DEFINE_string(galaxy_file_a, "",
              "Name of the ASCII file containing the first galaxies");
DEFINE_string(galaxy_file_b, "",
              "Name of the ASCII file containing the second galaxies");
DEFINE_bool(galaxy_radec, false, "Galaxy coordinates are in RA-DEC");
DEFINE_bool(use_only_pairs, false, "Use Only Pairs in correlation");
DEFINE_string(output_tag, "test",
              "Tag for output file: Wtheta_OUTPUT_TAG");
DEFINE_double(theta_min, 0.001, "Minimum angular scale (in degrees)");
DEFINE_double(theta_max, 10.0, "Maximum angular scale (in degrees)");
DEFINE_double(mag_min_a, 18.0, "Minimum acceptable galaxy magnitude");
DEFINE_double(mag_max_a, 28.0, "Maximum acceptable galaxy magnitude");
DEFINE_double(mag_min_b, 18.0, "Minimum acceptable galaxy magnitude");
DEFINE_double(mag_max_b, 28.0, "Maximum acceptable galaxy magnitude");
DEFINE_double(prob_min, 0.2, "Minimum acceptable galaxy likelihood");
DEFINE_double(prob_max, 1.00001, "Maximum acceptable galaxy likelihood");
DEFINE_int32(n_bins_per_decade, 5, "Number of angular bins per decade.");
DEFINE_int32(n_random, 1,
	     "Integer number of random points per galaxy to use.");
DEFINE_int32(n_jackknife, 0,
	     "Number of jack-knife samples to use. Defaults to 2*angular bins");
DEFINE_bool(single_index, false, "Use older single-index file format.");
DEFINE_bool(no_weight, false, "Map file is missing weight column.");
DEFINE_bool(coordinates_only, false, "Galaxy files only contain coordinates.");
DEFINE_int32(maximum_resolution, -1,
	     "Maximum resolution to use for pixel-based estimator");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --map_file=<StompMap ASCII>";
  usage += " --galaxy_file_a=<list of galaxy input files>";
  usage += " --galaxy_file_b=<list of galaxy input files>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_map_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_galaxy_file_a.empty() || FLAGS_galaxy_file_b.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  // First, we read our STOMP map into a map object.  There are a couple
  // permutations based on the various map formats that are out there: with
  // or without a weight column or in the single index or double index format.
  Stomp::Map* stomp_map;
  if (FLAGS_single_index) {
    if (FLAGS_no_weight) {
      stomp_map = new Stomp::Map(FLAGS_map_file, false, false);
    } else {
      stomp_map = new Stomp::Map(FLAGS_map_file, false);
    }
  } else {
    if (FLAGS_no_weight) {
      stomp_map = new Stomp::Map(FLAGS_map_file, true, false);
    } else {
      stomp_map = new Stomp::Map(FLAGS_map_file);
    }
  }
  std::cout << "Read map from " << FLAGS_map_file << "; total area: " <<
    stomp_map->Area() << " sq. deg.\n";


  // Now we read in our galaxy data file.  The expected format is
  //  LAMBDA  ETA  WEIGHT  MAGNITUDE
  // where the WEIGHT column is the likelihood that the object is a galaxy
  // and MAGNITUDE is the apparent magnitude in a given filter.  We filter all
  // of the objects against the map, tossing out any objects that aren't in the
  // map.
  Stomp::WAngularVector galaxy_a;
  Stomp::WAngularVector galaxy_b;

  std::cout << "Parsing " << FLAGS_galaxy_file_a << " and " << FLAGS_galaxy_file_b << " files...\n";
  unsigned long n_galaxy_a = 0;
  unsigned long n_galaxy_b = 0;

  Stomp::AngularCoordinate::Sphere galaxy_sphere =
    Stomp::AngularCoordinate::Survey;
  if (FLAGS_galaxy_radec) galaxy_sphere = Stomp::AngularCoordinate::Equatorial;
  std::ifstream galaxy_file_a(FLAGS_galaxy_file_a.c_str());
  std::ifstream galaxy_file_b(FLAGS_galaxy_file_b.c_str());
  double theta, phi, prob, mag, weight;
  
  prob = 1.0;
  mag = 0.5*(FLAGS_mag_max_a + FLAGS_mag_min_a);
  
  while (!galaxy_file_a.eof()) {
    char c;
    c = galaxy_file_a.peek();
    if (c == '#') {
      galaxy_file_a.ignore(2048, '\n');
      continue;
    }
    if (FLAGS_coordinates_only) {
      galaxy_file_a >> theta >> phi;
    } else {
      galaxy_file_a >> theta >> phi >> prob >> mag;
    }
    
    if (!galaxy_file_a.eof()) {
      Stomp::WeightedAngularCoordinate tmp_ang(theta, phi,
					       prob, galaxy_sphere);
      
      if ((mag >= FLAGS_mag_min_a) && (mag <= FLAGS_mag_max_a) &&
	  (stomp_map->FindLocation(tmp_ang,weight))) galaxy_a.push_back(tmp_ang);
      n_galaxy_a++;
    }
  }
  galaxy_file_a.close();
  std::cout << "Read " << n_galaxy_a << " galaxies; kept " <<
    galaxy_a.size() << "\n";
  n_galaxy_a = galaxy_a.size();
  galaxy_a.resize(n_galaxy_a);

  prob =1.0;
  mag = 0.5*(FLAGS_mag_max_b + FLAGS_mag_min_b);

  while (!galaxy_file_b.eof()) {
    char c;
    c = galaxy_file_b.peek();
    if (c == '#') {
      galaxy_file_b.ignore(2048, '\n');
      continue;
    }
    if (FLAGS_coordinates_only) {
      galaxy_file_b >> theta >> phi;
    } else {
      galaxy_file_b >> theta >> phi >> prob >> mag;
    }
    
    if (!galaxy_file_b.eof()) {
      Stomp::WeightedAngularCoordinate tmp_ang(theta, phi,
					       prob, galaxy_sphere);
      
      if ((prob >= FLAGS_prob_min) && (prob <= FLAGS_prob_max) &&
	  (mag >= FLAGS_mag_min_b) && (mag <= FLAGS_mag_max_b) &&
	  (stomp_map->FindLocation(tmp_ang,weight))) galaxy_b.push_back(tmp_ang);
      n_galaxy_b++;
    }
  }
  galaxy_file_b.close();
  std::cout << "Read " << n_galaxy_b << " galaxies; kept " <<
    galaxy_b.size() << "\n";
  n_galaxy_b = galaxy_b.size();
  galaxy_b.resize(n_galaxy_b);

  // Now, we set up the object that will contain the measurement results.  The
  // correlation object is a essentially a container for angular bin objects
  // which have a given angular range (all object or pixel pairs separated by
  // 0.01 < theta < 0.1 degrees, for instance).  In addition, the constructor
  // for these objects will work out, based on the angular bin size, which
  // StompMap resolution would be appropriate for calculating the angular
  // correlation on that scale.
  Stomp::AngularCorrelation wtheta(FLAGS_theta_min, FLAGS_theta_max,
				   FLAGS_n_bins_per_decade);
  // That pixel-based estimator works well on large scales, but on small scales
  // we want to use a pair-based estimator (which will be faster and require
  // less memory, provided we choose the break sensibly).  This call will
  // modify all of the high-resolution bins so that they use the pair-based
  // estimator.

  // Now we use the regions version of the auto-correlation code to find our
  // result.
  if (FLAGS_maximum_resolution == -1) {
    wtheta.AutoMaxResolution(static_cast<uint32_t>(sqrt(1.0*n_galaxy_a*
							n_galaxy_b)), 
			     stomp_map->Area());
  }
  else {
    std::cout << "Setting maximum resolution to " <<
      static_cast<uint16_t>(FLAGS_maximum_resolution) << "...\n";
    wtheta.SetMaxResolution(static_cast<uint16_t>(FLAGS_maximum_resolution));
  }
  if (FLAGS_use_only_pairs) {
    wtheta.UseOnlyPairs();
  }

  // Now we use the regions version of the auto-correlation code to find our
  // result.
  wtheta.FindCrossCorrelationWithRegions(*stomp_map, 
					galaxy_a,
					galaxy_b,
					static_cast<uint8_t>(FLAGS_n_random),
					static_cast<uint16_t>(FLAGS_n_jackknife
							      ));


  // And write out the results...
  std::string wtheta_file_name = "Wtheta_" + FLAGS_output_tag;
  std::cout << "Writing galaxy cross-correlation to " <<
    wtheta_file_name << "\n";
  wtheta.Write(wtheta_file_name);
  std::string wcov_file_name = "Wcov_" + FLAGS_output_tag;
  std::cout << "Writing galaxy cross-correlation covariance to " <<
    wcov_file_name << "\n";
  wtheta.WriteCovariance(wcov_file_name);

  return 0;
}
