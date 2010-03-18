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
DEFINE_string(background_files, "",
              "CSV names of the ASCII file containing the background objects");
DEFINE_bool(background_radec, false, "Background coordinates are in RA-DEC");
DEFINE_string(foreground_file, "",
              "ASCII file containing the locations to sample.");
DEFINE_bool(foreground_radec, false, "Foreground coordinates are in RA-DEC");
DEFINE_string(output_tag, "test",
              "Tag for output file: LocalDensity_OUTPUT_TAG");
DEFINE_double(r_min, 0.001, "Minimum projected scale (in Mpc/h)");
DEFINE_double(r_max, 100.0, "Maximum projected scale (in Mpc/h)");
DEFINE_int32(n_bins_per_decade, 4, "Number of angular bins per decade.");
DEFINE_bool(single_index, false, "Use older single-index file format.");
DEFINE_bool(no_weight, false, "Input file is missing weight column.");
DEFINE_int32(maximum_resolution, 128,
	     "Maximum resolution to use for pixel-based estimator");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --map_file=<StompMap ASCII>";
  usage += " --background_file=<Background catalog ASCII>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_map_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_background_files.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_foreground_file.empty()) {
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


  // Now we read in our background data files.  The expected format is
  //  LAMBDA  ETA  WEIGHT  MAGNITUDE
  // where the WEIGHT column is the likelihood that the object is a galaxy
  // and MAGNITUDE is the apparent magnitude in a given filter.  We filter all
  // of the objects against the map, tossing out any objects that aren't in the
  // map.
  Stomp::WAngularVector background;

  // First we extract the background file names into a vector of strings.  The
  // input vector of file names should be comma-separated.
  std::vector<std::string> background_files;

  Stomp::Tokenize(FLAGS_background_files, background_files, ",");

  std::cout << "Parsing " << background_files.size() << " files...\n";
  unsigned long n_background = 0;

  // Now iterate through our list of files, appending the galaxies if they are
  // in the input map.
  Stomp::AngularCoordinate::Sphere background_sphere =
    Stomp::AngularCoordinate::Survey;
  if (FLAGS_background_radec)
    background_sphere = Stomp::AngularCoordinate::Equatorial;
  for (std::vector<std::string>::iterator file_iter=background_files.begin();
       file_iter!=background_files.end();++file_iter) {

    std::cout << "\tParsing " << file_iter->c_str() << "...\n";
    std::ifstream background_file(file_iter->c_str());
    double lambda, eta, prob, mag;

    while (!background_file.eof()) {
      background_file >> lambda >> eta >> prob >> mag;
      Stomp::WeightedAngularCoordinate tmp_ang(lambda, eta, prob,
					       background_sphere);

      if (stomp_map->FindLocation(tmp_ang) &&
	  (tmp_ang.Weight() > 0.2)) background.push_back(tmp_ang);
      n_background++;
    }
    background_file.close();
  }

  std::cout << "Read " << n_background <<
    " galaxies from input files; kept " << background.size() << "\n\n";
  n_background = background.size();

  // Now we read in the foreground objects.
  Stomp::WAngularVector foreground;

  std::cout << "Parsing " << FLAGS_foreground_file << "...\n";
  unsigned long n_foreground = 0;

  Stomp::AngularCoordinate::Sphere foreground_sphere =
    Stomp::AngularCoordinate::Survey;
  if (FLAGS_foreground_radec)
    foreground_sphere = Stomp::AngularCoordinate::Equatorial;
  std::cout << "\tParsing " << FLAGS_foreground_file << "...\n";
  std::ifstream foreground_file(FLAGS_foreground_file.c_str());
  double theta, phi, weight;

  while (!foreground_file.eof()) {
    foreground_file >> theta >> phi >> weight;
    Stomp::WeightedAngularCoordinate tmp_ang(theta, phi, weight,
					     foreground_sphere);

    if (stomp_map->FindLocation(tmp_ang)) foreground.push_back(tmp_ang);
    n_foreground++;
  }

  foreground_file.close();

  std::cout << "Read " << n_foreground <<
    " foreground objects from input file; kept " << foreground.size() << "\n";
  n_foreground = foreground.size();

  // Normally, the AngularCorrelation object is used to store angular-binned
  // measurements, but we can use it to set up our fidducial projected radial
  // binning.  The minimum and maximum radii are taken to be in Mpc/h.
  Stomp::AngularCorrelation wrp(FLAGS_r_min, FLAGS_r_max,
				FLAGS_n_bins_per_decade);

  // We want to find the mean weight of our background sample in the same
  // radial bins around each foreground object.  This is complicated by the
  // fact that the foreground objects are at different redshifts.
  //
  // Before getting to that, though, we'll want to create a TreeMap and
  // populate it with our background objects.  The resolution of the tree
  // is somewhat arbitrary.
  std::cout << "\tBuilding background tree...\n";
  Stomp::TreeMap* background_tree =
    new Stomp::TreeMap(wrp.MinResolution(), 200);

  for (Stomp::WAngularIterator iter=background.begin();
       iter!=background.end();++iter) {
    if (!background_tree->AddPoint(*iter)) {
      std::cout << "Failed to add point: " << iter->Lambda() << ", " <<
	iter->Eta() << "\n";
    }
  }

  // Now we can clear the background objects from memory
  background.clear();

  // Now we iterate over the foregrounds and angular bins...
  for (Stomp::WAngularIterator foreground_iter=foreground.begin();
       foreground_iter!=foreground.end();++foreground_iter) {
    // For a given foreground object, we need to set up an AngularCorrelation
    // object based on the redshift of the foreground object, which we'll take
    // to have been stored in that object's Weight() value.
    Stomp::AngularCorrelation tmp_wtheta = wrp;

    // Now we need to iterate through the radial bins and find their angular
    // extents.
    for (unsigned long i=0;i<wrp.NBins();i++) {
      Stomp::ThetaIterator rp_iter = wrp.BinIterator(i);
      Stomp::ThetaIterator tmp_iter = tmp_wtheta.BinIterator(i);

      // Now we translate the radial bins to angular scales.
      tmp_iter->
	SetThetaMin(Stomp::Cosmology::ProjectedAngle(foreground_iter->Weight(),
						     rp_iter->ThetaMin()));
      tmp_iter->
	SetThetaMax(Stomp::Cosmology::ProjectedAngle(foreground_iter->Weight(),
						     rp_iter->ThetaMax()));
    }

    // Now we iterate over our temporary AngularCorrelation object, finding
    // pairs.
    for (Stomp::ThetaIterator theta_iter=tmp_wtheta.Begin();
	 theta_iter!=tmp_wtheta.End();++theta_iter) {
      theta_iter->Reset();

      // This will store the sum of the weights in the Weight field and
      // the number of pairs in the Counter field for each angular bin.
      background_tree->FindWeightedPairs(*foreground_iter, *theta_iter);
    }

    // Once we've iterated through all of the angular bins, we add the
    // values from each bin to our radial bin version
    for (unsigned long i=0;i<wrp.NBins();i++) {
      Stomp::ThetaIterator rp_iter = wtheta.BinIterator(i);
      Stomp::ThetaIterator tmp_iter = tmp_wtheta.BinIterator(i);

      rp_iter->AddToWeight(tmp_iter->Weight());
      rp_iter->AddToCounter(tmp_iter->Counter());
    }
  }

  // Finally, we write out the results...
  std::string proj_file_name = "ProjectedDensity_" + FLAGS_output_tag;
  std::cout << "Writing projected densities to " <<
    proj_file_name << "\n";

  std::ofstream output_file(proj_file_name.c_str());
  for (int i=0;i<foreground.size();i++) {
    // First write out the foreground weight.
    output_file << foreground[i].Weight();

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
