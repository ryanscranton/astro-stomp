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
DEFINE_string(galaxy_files, "",
              "CSV names of the ASCII file containing the galaxy catalog");
DEFINE_bool(galaxy_radec, false, "Galaxy coordinates are in RA-DEC");
DEFINE_string(target_file, "",
              "ASCII file containing the locations to sample.");
DEFINE_bool(target_radec, false, "Target coordinates are in RA-DEC");
DEFINE_string(output_tag, "test",
              "Tag for output file: LocalDensity_OUTPUT_TAG");
DEFINE_double(theta_min, 0.001, "Minimum angular scale (in degrees)");
DEFINE_double(theta_max, 1.0, "Maximum angular scale (in degrees)");
DEFINE_int32(n_bins_per_decade, 4, "Number of angular bins per decade.");
DEFINE_bool(single_index, false, "Use older single-index file format.");
DEFINE_bool(no_weight, false, "Input file is missing weight column.");
DEFINE_int32(maximum_resolution, 128,
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

  if (FLAGS_galaxy_files.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_target_file.empty()) {
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


  // Now we read in our galaxy data files.  The expected format is
  //  LAMBDA  ETA  WEIGHT  MAGNITUDE
  // where the WEIGHT column is the likelihood that the object is a galaxy
  // and MAGNITUDE is the apparent magnitude in a given filter.  We filter all
  // of the objects against the map, tossing out any objects that aren't in the
  // map.
  Stomp::WAngularVector galaxy;

  // First we extract the galaxy file names into a vector of strings.  The
  // input vector of file names should be comma-separated.
  std::vector<std::string> galaxy_files;

  Stomp::Tokenize(FLAGS_galaxy_files, galaxy_files, ",");

  std::cout << "Parsing " << galaxy_files.size() << " files...\n";
  unsigned long n_galaxy = 0;

  // Now iterate through our list of files, appending the galaxies if they are
  // in the input map.
  Stomp::AngularCoordinate::Sphere galaxy_sphere =
    Stomp::AngularCoordinate::Survey;
  if (FLAGS_galaxy_radec) galaxy_sphere = Stomp::AngularCoordinate::Equatorial;
  for (std::vector<std::string>::iterator file_iter=galaxy_files.begin();
       file_iter!=galaxy_files.end();++file_iter) {

    std::cout << "\tParsing " << file_iter->c_str() << "...\n";
    std::ifstream galaxy_file(file_iter->c_str());
    double lambda, eta, prob, mag;

    while (!galaxy_file.eof()) {
      galaxy_file >> lambda >> eta >> prob >> mag;
      Stomp::WeightedAngularCoordinate tmp_ang(lambda, eta, prob,
					       galaxy_sphere);

      if (stomp_map->FindLocation(tmp_ang) &&
	  (tmp_ang.Weight() > 0.2)) galaxy.push_back(tmp_ang);
      n_galaxy++;
    }
    galaxy_file.close();
  }

  std::cout << "Read " << n_galaxy << " galaxies from input files; kept " <<
    galaxy.size() << "\n\n";
  n_galaxy = galaxy.size();

  // Now we read in the target objects.
  Stomp::WAngularVector target;

  std::cout << "Parsing " << FLAGS_target_file << "...\n";
  unsigned long n_target = 0;

  Stomp::AngularCoordinate::Sphere target_sphere =
    Stomp::AngularCoordinate::Survey;
  if (FLAGS_target_radec) target_sphere = Stomp::AngularCoordinate::Equatorial;
  std::cout << "\tParsing " << FLAGS_target_file << "...\n";
  std::ifstream target_file(FLAGS_target_file.c_str());
  double theta, phi, weight;

  while (!target_file.eof()) {
    target_file >> theta >> phi >> weight;
    Stomp::WeightedAngularCoordinate tmp_ang(theta, phi, weight,
					     target_sphere);

    if (stomp_map->FindLocation(tmp_ang)) target.push_back(tmp_ang);
    n_target++;
  }

  target_file.close();

  std::cout << "Read " << n_target << " targets from input file; kept " <<
    target.size() << "\n";
  n_target = target.size();

  // Now, we set up the object that will contain the measurement results.  The
  // correlation object is a essentially a container for angular bin objects
  // which have a given angular range (all object or pixel pairs separated by
  // 0.01 < theta < 0.1 degrees, for instance).  In addition, the constructor
  // for these objects will work out, based on the angular bin size, which
  // Stomp::Map resolution would be appropriate for calculating the angular
  // correlation on that scale.
  Stomp::AngularCorrelation wtheta(FLAGS_theta_min, FLAGS_theta_max,
				   FLAGS_n_bins_per_decade);

  // We want to find the local galaxy density in annuli around each of our
  // target points.  To do this, we need to know how many galaxies are in each
  // of those annuli as well as the area of each annulus within the survey
  // area.
  //
  // For the first task, we'll want to create a TreeMap and populate
  // it with our galaxies.  The resolution of the tree is somewhat arbitrary.
  std::cout << "\tBuilding galaxy tree...\n";
  Stomp::TreeMap* galaxy_tree = new Stomp::TreeMap(wtheta.MinResolution(), 200);

  for (Stomp::WAngularIterator iter=galaxy.begin();iter!=galaxy.end();++iter) {
    if (!galaxy_tree->AddPoint(*iter)) {
      std::cout << "Failed to add point: " << iter->Lambda() << ", " <<
	iter->Eta() << "\n";
    }
  }

  // Now we can clear the galaxies from memory
  galaxy.clear();

  // We're going to want to store a copy of the AngularCorrelation object for
  // each of the targets, so we'll set up that vector now.
  Stomp::WThetaVector wthetaVec;
  wthetaVec.reserve(target.size());

  // Now we iterate over the targets and angular bins...
  for (Stomp::WAngularIterator target_iter=target.begin();
       target_iter!=target.end();++target_iter) {
    for (Stomp::ThetaIterator theta_iter=wtheta.Begin();
	 theta_iter!=wtheta.End();++theta_iter) {
      theta_iter->Reset();

      // For each angular bin, we find the pixels at the corresponding
      // resolution that cover the bin.
      Stomp::Pixel center_pix(*target_iter, theta_iter->Resolution());
      Stomp::PixelVector pixVec;
      center_pix.WithinAnnulus(*theta_iter, pixVec);

      // For each of those pixels, we'll want to find the fraction of the pixel
      // that's within our survey bounds as well as the number of galaxies that
      // are in the pixel.
      double total_area = 0.0, total_weight = 0.0;
      for (Stomp::PixelIterator pix_iter=pixVec.begin();
	   pix_iter!=pixVec.end();++pix_iter) {
	total_area += stomp_map->FindUnmaskedFraction(*pix_iter);
	total_weight += galaxy_tree->Weight(*pix_iter);
      }

      // Now we store that value in the corresponding angular bin's Weight
      // field.
      theta_iter->AddToWeight(total_weight/(total_area*center_pix.Area()));
    }

    // Once we've iterated through all of the angular bins, we append a copy
    // of the AngularCorrelation to our output vector to store the results for
    // this target.
    wthetaVec.push_back(wtheta);
  }

  // Finally, we write out the results...
  std::string wtheta_file_name = "LocalDensity_" + FLAGS_output_tag;
  std::cout << "Writing local densities to " <<
    wtheta_file_name << "\n";

  std::ofstream output_file(wtheta_file_name.c_str());
  for (int i=0;i<target.size();i++) {
    // First write out the target weight.
    output_file << target[i].Weight();

    // Now the densities for each angular bin...
    for (Stomp::ThetaIterator iter=wthetaVec[i].Begin();
	 iter!=wthetaVec[i].End();++iter) {
      output_file << " " << std::setprecision(6) << iter->Weight();
    }
    output_file << "\n";
  }
  output_file.close();

  return 0;
}
