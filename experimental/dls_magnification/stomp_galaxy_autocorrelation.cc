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
DEFINE_string(galaxy_file, "",
              "Name of the ASCII file containing the input galaxy catalog");
DEFINE_string(output_tag, "test",
              "Tag for output file: Wtheta_OUTPUT_TAG");
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

  if (FLAGS_galaxy_file.empty()) {
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
  //  RA  DEC  WEIGHT  MAGNITUDE
  // where the WEIGHT column is the likelihood that the object is a galaxy
  // and MAGNITUDE is the apparent magnitude in a given filter.  We filter all
  // of the objects against the map, tossing out any objects that aren't in the
  // map.
  Stomp::WAngularVector galaxy;
  std::ifstream galaxy_file(FLAGS_galaxy_file.c_str());
  double ra, dec, prob, mag;
  uint32_t n_galaxy = 0;
  
  std::cout << "Reading in File\n";
  while (!galaxy_file.eof()) {
    galaxy_file >> ra >> dec >> prob >> mag;
    Stomp::WeightedAngularCoordinate tmp_ang(ra, dec, prob,
					     Stomp::AngularCoordinate::Equatorial);
    if (stomp_map->FindLocation(tmp_ang) &&
	(tmp_ang.Weight() > 0.2)) galaxy.push_back(tmp_ang);
    n_galaxy++;
  }
  galaxy_file.close();

  std::cout << "Read " << n_galaxy << " galaxies from " << FLAGS_galaxy_file <<
    "; kept " << galaxy.size() <<"\n";
  n_galaxy = galaxy.size();


  // Now, we set up the object that will contain the measurement results.  The
  // correlation object is a essentially a container for angular bin objects
  // which have a given angular range (all object or pixel pairs separated by
  // 0.01 < theta < 0.1 degrees, for instance).  In addition, the constructor
  // for these objects will work out, based on the angular bin size, which
  // Stomp::Map resolution would be appropriate for calculating the angular
  // correlation on that scale.
  Stomp::AngularCorrelation wtheta(FLAGS_theta_min, FLAGS_theta_max,
				   FLAGS_n_bins_per_decade);

  // That pixel-based estimator works well on large scales, but on small scales
  // we want to use a pair-based estimator (which will be faster and require
  // less memory, provided we choose the break sensibly).  This call will
  // modify all of the high-resolution bins so that they use the pair-based
  // estimator.
  wtheta.SetMaxResolution(FLAGS_maximum_resolution);


  // Now, we're ready to start calculating the autocorrelation.  It is
  // possible to do this all in a single call with
  //
  // wtheta.FindAutoCorrelation(*stomp_map, galaxy, 1);
  //
  // where the last argument controls the number of random catalogs generated
  // for the pair-based estimator.  Instead, we'll walk through the process
  // explicitly, starting with the creation of a ScalarMap for the
  // pixel-based estimator.
  std::cout << "Using pixel-based estimator for " <<
    wtheta.ThetaMin(wtheta.MaxResolution()) << " < theta < " <<
    wtheta.ThetaMax(wtheta.MinResolution()) << "...\n";

  std::cout << "\tMaking density map...\n";
  Stomp::ScalarMap* scalar_map =
    new Stomp::ScalarMap(*stomp_map, FLAGS_maximum_resolution,
			 Stomp::ScalarMap::DensityField);
  n_galaxy = 0;
  for (Stomp::WAngularIterator iter=galaxy.begin();iter!=galaxy.end();++iter)
    if (scalar_map->AddToMap(*iter)) n_galaxy++;

  // From here, we could use the FindPixelAutoCorrelation() method on the
  // AngularCorrelation object:
  //
  // wtheta.FindPixelAutoCorrelation(*density_map);
  //
  // but we'll do this by hand for demonstration purposes.
  std::cout << "\tAuto-correlating...\n";
  scalar_map->AutoCorrelate(wtheta);

  // Now, we iterate through our resolutions, making re-sampled maps and
  // calculating auto-correlations for the corresponding bins.
  for (uint16_t resolution=scalar_map->Resolution()/2;
       resolution>=wtheta.MinResolution();resolution/=2) {
    Stomp::ScalarMap* sub_scalar_map =
      new Stomp::ScalarMap(*scalar_map, resolution);
    sub_scalar_map->AutoCorrelate(wtheta);

    delete sub_scalar_map;
  }
  delete scalar_map;


  // We can dump our density map now and start doing the pair-based estimator.
  // By using the "0" value as the argument to the Theta method, we
  // are selecting the bins where the pair-based estimator has been used.
  std::cout << "Using pair-based estimator for " <<
    wtheta.ThetaMin(0) << " < theta < " << wtheta.ThetaMax(0) << "...\n";

  // As with the pixel-based estimator, this could all be done with a single
  // command:
  //
  // wtheta.FindPairAutoCorrelation(*stomp_map, galaxy);
  //
  // but we'll do things the long way to illustrate what's going on.
  std::cout << "\tBuilding galaxy tree...\n";
  Stomp::TreeMap* galaxy_tree =
    new Stomp::TreeMap(wtheta.MaxResolution(), 200);

  for (Stomp::WAngularIterator iter=galaxy.begin();iter!=galaxy.end();++iter) {
    if (!galaxy_tree->AddPoint(*iter)) {
      std::cout << "Failed to add point: " << iter->Lambda() << ", " <<
	iter->Eta() << "\n";
    }
  }


  // Galaxy-galaxy
  std::cout << "\tCalculating galaxy-galaxy pairs...\n";
  galaxy_tree->FindWeightedPairs(galaxy, wtheta);

  // After finding the weighted pairs, we need to transfer the results from
  // the weight field in each angular bin to the galaxy-galaxy field.
  for (Stomp::ThetaIterator iter=wtheta.Begin(0);iter!=wtheta.End(0);++iter)
    iter->MoveWeightToGalGal();

  std::cout << "\tGenerating random points...\n";
  Stomp::WAngularVector random_galaxy;
  stomp_map->GenerateRandomPoints(random_galaxy, galaxy);
  // We can clear the galaxy data out of memory now.
  galaxy.clear();

  // Galaxy-Random -- there's a symmetry here, so the results go in GalRand
  // and RandGal.
  std::cout << "\tCalculating galaxy-random pairs...\n";
  galaxy_tree->FindWeightedPairs(random_galaxy, wtheta);
  for (Stomp::ThetaIterator iter=wtheta.Begin(0);iter!=wtheta.End(0);++iter)
    iter->MoveWeightToGalRand(true);

  // We're done with the galaxy tree, so we can dump it from memory
  delete galaxy_tree;

  // ... and build the random point tree.
  std::cout << "\tBuilding random galaxy tree...\n";
  Stomp::TreeMap* random_tree =
    new Stomp::TreeMap(Stomp::HPixResolution, 200);

  for (Stomp::WAngularIterator iter=random_galaxy.begin();
       iter!=random_galaxy.end();++iter) {
    if (!random_tree->AddPoint(*iter)) {
      std::cout << "Failed to add point: " << iter->Lambda() << ", " <<
	iter->Eta() << "\n";
    }
  }

  // Random-Random
  std::cout << "\tCalculating random-random pairs...\n";
  random_tree->FindWeightedPairs(random_galaxy, wtheta);
  for (Stomp::ThetaIterator iter=wtheta.Begin(0);iter!=wtheta.End(0);++iter)
    iter->MoveWeightToRandRand();

  // Finally write out the results...
  std::string wtheta_file_name = "Wtheta_" + FLAGS_output_tag;
  std::cout << "Writing galaxy auto-correlation to " <<
    wtheta_file_name << "\n";

  std::ofstream output_file(wtheta_file_name.c_str());
  for (Stomp::ThetaIterator iter=wtheta.Begin();iter!=wtheta.End();++iter) {
    output_file << std::setprecision(6) << iter->Theta() << " " <<
      iter->Wtheta()  << " " << iter->WthetaError() << "\n";
  }
  output_file.close();

  return 0;
}
