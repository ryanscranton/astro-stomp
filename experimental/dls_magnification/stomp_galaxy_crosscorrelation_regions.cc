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
DEFINE_string(region_file, "",
	      "ASCII file containing region definitions");
DEFINE_bool(galaxy_radec, false, "Galaxy coordinates are in RA-DEC");
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
DEFINE_int32(n_jack, -1, 
	     "Number of Jackknife Samples. Defaults to number of regions specifitied by the region file if not set or less than region file.");
DEFINE_int32(n_random, 1,
	     "Integer number of random points per galaxy to use.");
DEFINE_bool(single_index, false, "Use older single-index file format.");
DEFINE_bool(no_weight, false, "Map file is missing weight column.");
DEFINE_bool(coordinates_only, false, "Galaxy files only contain coordinates.");
DEFINE_int32(maximum_resolution, -1,
	     "Maximum resolution to use for pixel-based estimator");
DEFINE_string(density_output,"testDensity.map",
	      "Get Rid of this later");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --map_file=<StompMap ASCII>";
  usage += " --galaxy_file_a=<list of galaxy input files>";
  usage += " --galaxy_file_b=<list of galaxy input files>";
  usage += " --region_file=<ASCII region definitions>";
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
  if (FLAGS_region_file.empty()) {
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

  int type;
  int n_regions = 0;
  double region_area =0;
  double radius;
  double ra[4], dec[4];
  unsigned long idx;
  double weight;
  std::ifstream region_file(FLAGS_region_file.c_str());
  Stomp::RegionBoundVector region_vector;
  
  std::cout << "Reading regions from " << FLAGS_region_file; 
  while (!region_file.eof()) {
    region_file >> idx >> type >> radius >> weight >> ra[0] >> dec[0] >> ra[1] >> dec[1] >> ra[2] >> dec[2] >> ra[3] >> dec[3];
    
    if (!region_file.eof()) {
      
      Stomp::AngularVector ang_vector;
      for (int i=0;i<4;i++)
	ang_vector.push_back(Stomp::AngularCoordinate(ra[i],dec[i],Stomp::AngularCoordinate::Equatorial));
      
      Stomp::PolygonBound* poly_bound = new Stomp::PolygonBound(ang_vector);
      Stomp::RegionBound region_bound(poly_bound);
      region_area += region_bound.BoundArea();
      
      region_vector.push_back(region_bound);
      n_regions++;
    }
  }
  std::cout << "; " << n_regions << " Loaded; Total Area " << region_area 
	    << "; Average Area " << region_area/(1.0*n_regions) << "\n";

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
  double theta, phi, prob, mag;
  
  prob = 1.0;
  mag = 0.5*(FLAGS_mag_max_a + FLAGS_mag_min_a);
  
  while (!galaxy_file_a.eof()) {
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

  mag = 0.5*(FLAGS_mag_max_b + FLAGS_mag_min_b);

  while (!galaxy_file_b.eof()) {
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
  wtheta.AutoMaxResolution(sqrt(1.0*n_galaxy_a*n_galaxy_b),stomp_map->Area());
  //else {
  //  std::cout << "Setting maximum resolution to " <<
  //    static_cast<int>(FLAGS_maximum_resolution) << "...\n";
  //  wtheta.SetMaxResolution(static_cast<uint16_t>(FLAGS_maximum_resolution));
  //}
  std::cout << "Maximum resolution Set to " <<
    static_cast<int>(wtheta.MaxResolution()) << "...\n";

  // That pixel-based estimator works well on large scales, but on small scales
  // we want to use a pair-based estimator (which will be faster and require
  // less memory, provided we choose the break sensibly).  This call will
  // modify all of the high-resolution bins so that they use the pair-based
  // estimator.
  std::cout << "Using pixel-based estimator for " <<
    wtheta.ThetaMin(wtheta.MaxResolution()) << " < theta < " <<
    wtheta.ThetaMax(wtheta.MinResolution()) << "...\n";

  Stomp::PixelVector coverage_pix;
  Stomp::PixelVector hold_pix;
  stomp_map->Coverage(coverage_pix, wtheta.MinResolution());
  for (Stomp::PixelIterator iter=coverage_pix.begin();
       iter!=coverage_pix.end();++iter) {
    for (uint32_t bound_iter=0;bound_iter<region_vector.size();bound_iter++) {
      if (region_vector[bound_iter].CheckPixel(*iter)) {
	hold_pix.push_back(*iter);
	break;
      }
    }
  }
  Stomp::Map* temp_map = new Stomp::Map(hold_pix);
  stomp_map->IntersectMap(*temp_map);

  std::cout<< "Initializing Regions....\n";
  if (FLAGS_n_jack>-1 && FLAGS_n_jack>n_regions)
    n_regions = FLAGS_n_jack;
  uint16_t n_regions_final =
    stomp_map->InitializeRegions(region_vector, n_regions, wtheta.MinResolution());
  std::cout << "Requested " << n_regions << " regions; ended up with " << n_regions_final << "...\n";

  for (Stomp::ThetaIterator iter=wtheta.Begin();iter!=wtheta.End();++iter)
    iter->InitializeRegions(n_regions_final);

  std::cout<< "Generating Coverage Map...\n";
  //Stomp::PixelVector coverage_map;
  //stomp_map->Coverage(coverage_map, wtheta.MaxResolution());
  Stomp::Map* density_map = 
    new Stomp::Map(*stomp_map);
  density_map->InitializeRegions(*stomp_map);

  std::cout<< "Creating Scalar Maps at " << wtheta.MaxResolution() << "...\n";
  Stomp::ScalarMap* galaxy_map_a =
    new Stomp::ScalarMap(*stomp_map, wtheta.MaxResolution(),
			 Stomp::ScalarMap::DensityField);
  
  Stomp::ScalarMap* galaxy_map_b =
    new Stomp::ScalarMap(*galaxy_map_a, wtheta.MaxResolution(),
			 Stomp::ScalarMap::DensityField);
  std::cout<< "Populating Scalar Maps...\n";
  n_galaxy_a = 0;
  for (Stomp::WAngularIterator iter=galaxy_a.begin();iter!=galaxy_a.end();++iter)
    if (galaxy_map_a->AddToMap(*iter)) n_galaxy_a++;
  galaxy_map_a->InitializeRegions(*stomp_map);
  galaxy_map_a->UseLocalMeanIntensity(true);

  n_galaxy_b=0;
  for (Stomp::WAngularIterator iter=galaxy_b.begin();iter!=galaxy_b.end();++iter)
    if (galaxy_map_b->AddToMap(*iter)) n_galaxy_b++;
  galaxy_map_b->InitializeRegions(*stomp_map);
  galaxy_map_b->UseLocalMeanIntensity(true);

  std::cout << "Imprinting Stomp map with Scalar density.\n" 
	    << "NPoints is " << galaxy_map_b->NPoints() 
	    << "; Density is " << galaxy_map_b->Density() 
	    << "; Mean Intesity is " << galaxy_map_b->MeanIntensity() << "...\n";

  std::cout << "Creating Density Map...\n";
  galaxy_map_b->ImprintMap(*density_map, true);
  std::cout << "Density Map is " << density_map->Area() << " sq. deg; " << density_map->Size() << " Size; " << density_map->MinWeight() << " MinWeight; " << density_map->MaxWeight() << " MaxWeight;\n";

  double scale = 1.0/galaxy_map_b->MeanIntensity();
  double min_weight = density_map->MinWeight()*scale;
  double max_weight = density_map->MaxWeight()*scale;
  
  std::cout << "\tScaling Weight By Value " << scale << "; New Min and Max Should be " << min_weight << ", " << max_weight << "\n";
  density_map->ScaleWeight(scale);
  // tmp_map->Write(FLAGS_density_output);
  // Stomp::Map* density_map;
  // density_map = new Stomp::Map(FLAGS_density_output);
  // density_map->InitializeRegions(*stomp_map);
  //Stomp::Map* density_map = 
  //  new Stomp::Map(*tmp_map);
  std::cout << "Density Map is " << density_map->Area() << " sq. deg; " << density_map->Size() << " Size; " << density_map->MinWeight() << " MinWeight; " << density_map->MaxWeight() << " MaxWeight;\n";
  // density_map->InitializeRegions(*stomp_map);
  // density_map->Write(FLAGS_density_output);
  
  std::cout << "\tCross-correlating...\n";
  galaxy_map_a->CrossCorrelateWithRegions(*galaxy_map_b, wtheta);
  for (uint32_t resolution=galaxy_map_a->Resolution()/2;
       resolution>=wtheta.MinResolution();resolution/=2) {
    
    std::cout << "\t\tworking on resolution "<< resolution <<"...\n";
    Stomp::ScalarMap* sub_galaxy_map_a =
      new Stomp::ScalarMap(*galaxy_map_a, resolution);
    sub_galaxy_map_a->InitializeRegions(*stomp_map);
    sub_galaxy_map_a->UseLocalMeanIntensity(true);
    
    Stomp::ScalarMap* sub_galaxy_map_b =
      new Stomp::ScalarMap(*galaxy_map_b, resolution);
    sub_galaxy_map_b->InitializeRegions(*stomp_map);
    sub_galaxy_map_b->UseLocalMeanIntensity(true);
    sub_galaxy_map_a->CrossCorrelateWithRegions(*sub_galaxy_map_b, wtheta);
    delete sub_galaxy_map_a;
    delete sub_galaxy_map_b;
  }
  std::cout << "Done with Pixel Based correlator...\n";
  delete galaxy_map_a;
  delete galaxy_map_b;

  //Now It's time to run through the pair based estimator to do this we need
  //to create the tree map and regionate as above but with the added twist of
  //using the previous scalar map to set up the local galaxy densities on each
  //region for use in the random tree map. right now this assumes that galaxy_map_b
  //is the higher density foreground map. Later interations of this code may have a
  //better way to compute the local denisty.

  //create our density map to produce randoms on


  //this should be mostly straight out of stomp_angular_correlation.cc
  std::cout << "Using pair-based estimator for " <<
    wtheta.ThetaMin(0) << " < theta < " << wtheta.ThetaMax(0) << "...\n";

  std::cout << "\tBuilding galaxy tree...\n";
  Stomp::TreeMap* galaxy_tree_a = new Stomp::TreeMap(wtheta.MinResolution(), 200);

  unsigned long int n_kept = 0;
  unsigned long int n_fail = 0;
  //load galaxies into tree
  for (Stomp::WAngularIterator iter=galaxy_a.begin();iter!=galaxy_a.end();++iter) {
    if (stomp_map->Contains(*iter)) {
      n_kept++;
      if (!galaxy_tree_a->AddPoint(*iter)) {
	std::cout << "Failed to add point: " << iter->Lambda() << ", " <<
	  iter->Eta() << "\n";
	n_fail++;
      }
    }
  }
  std::cout << n_kept - n_fail << "/" << galaxy_a.size() <<
    " objects added to tree;" << n_fail << " failed adds...\n";
  std::cout << "Nodes in map " << galaxy_tree_a->Nodes() << "; Resolution " <<
    galaxy_tree_a->Resolution() << "; NPoints" << galaxy_tree_a->NPoints() << "\n";

  
  //regionate galaxy_tree_a
  if (!galaxy_tree_a->InitializeRegions(*stomp_map)) {
    std::cout << "Failed to initialize regions on TreeMap  Exiting.\n";
    exit(2);
  } else {
    std::cout << "Regionated TreeMap with " << galaxy_tree_a->NRegion() << " regions...\n";
  }

  //Galaxy-Galaxy
  std::cout << "Working on Galaxy-Galaxy Pairs...\n";
  galaxy_tree_a->FindWeightedPairsWithRegions(galaxy_b, wtheta);
  // After finding the weighted pairs, we need to transfer the results from
  // the weight field in each angular bin to the galaxy-galaxy field.
  for (Stomp::ThetaIterator iter=wtheta.Begin(0);iter!=wtheta.End(0);++iter) {
    iter->MoveWeightToGalGal();
    std::cout << "Mean gal-gal is " << iter->MeanGalGal() << " " << iter->NRegion() << "\n";
  }

  // Before we start on the random iterations, we'll zero out the data fields
  // for those counts.
  //Time to generate random points and get this party started
  for (Stomp::ThetaIterator iter=wtheta.Begin(0);iter!=wtheta.End(0);++iter) {
    iter->ResetGalRand();
    iter->ResetRandGal();
    iter->ResetRandRand();
  }

  for (uint8_t rand_iter=0;rand_iter<FLAGS_n_random;rand_iter++) {
    std::cout << "\tRandom iteration " <<
      static_cast<int>(rand_iter) << "...\n";
    
    std::cout << "Generating random galaxy a...\n";
    Stomp::WAngularVector random_galaxy_a;
    density_map->GenerateRandomPoints(random_galaxy_a, galaxy_a, true);

    std::cout << "Generating random galaxy b...\n";
    Stomp::WAngularVector random_galaxy_b;
    density_map->GenerateRandomPoints(random_galaxy_b, galaxy_b, true);

    std::cout << "Length of random Galaxy Vectors " << random_galaxy_a.size() << "; " << random_galaxy_b.size() << "\n";

    //Galaxy_Random
    galaxy_tree_a->FindWeightedPairsWithRegions(random_galaxy_b, wtheta);
    for (Stomp::ThetaIterator iter=wtheta.Begin(0);iter!=wtheta.End(0);++iter) {
      //galaxy_tree_a->FindWeightedPairsWithRegions(random_galaxy_b, *iter);
      iter->MoveWeightToGalRand();
      std::cout << "Mean gal-rand is " << iter->MeanGalRand() << " " << iter->NRegion() << "\n";
    }

    Stomp::TreeMap* random_tree_a = new Stomp::TreeMap(wtheta.MaxResolution(), 200);

    for (Stomp::WAngularIterator iter=random_galaxy_a.begin();
	 iter!=random_galaxy_a.end();++iter) {
      if (!random_tree_a->AddPoint(*iter)) {
	std::cout << "Failed to add point: " << iter->Lambda() << ", " <<
	  iter->Eta() << "\n";
      }
    }
    if (!random_tree_a->InitializeRegions(*stomp_map)) {
      std::cout << "Failed to initialize regions on TreeMap  Exiting.\n";
      exit(2);
    }

    //Random-Galaxy
    random_tree_a->FindWeightedPairsWithRegions(galaxy_b, wtheta);
    for (Stomp::ThetaIterator iter=wtheta.Begin(0);iter!=wtheta.End(0);++iter) {
      //random_tree_a->FindWeightedPairsWithRegions(galaxy_b, *iter);
      iter->MoveWeightToRandGal();
      std::cout << "Mean rand-gal is " << iter->MeanRandGal() << " " << iter->NRegion() << "\n";
    }

    // Random-Random
    random_tree_a->FindWeightedPairsWithRegions(random_galaxy_b, wtheta);
    for (Stomp::ThetaIterator iter=wtheta.Begin(0);iter!=wtheta.End(0);++iter) {
      //random_tree_a->FindWeightedPairsWithRegions(random_galaxy_b, *iter);
      iter->MoveWeightToRandRand();
      std::cout << "Mean rand-rand is " << iter->MeanRandRand() << " " << iter->NRegion() << "\n";
    }

    delete random_tree_a;
  }

  delete galaxy_tree_a;

  // Finally, we rescale our random pair counts to normalize them to the
  // number of input objects.
  for (Stomp::ThetaIterator iter=wtheta.Begin(0);iter!=wtheta.End(0);++iter) {
    iter->RescaleGalRand(1.0*FLAGS_n_random);
    iter->RescaleRandGal(1.0*FLAGS_n_random);
    iter->RescaleRandRand(1.0*FLAGS_n_random);
  }

  // And write out the results...
  std::string wtheta_file_name = "Wtheta_" + FLAGS_output_tag;
  std::string wcovar_file_name = "Wcovar_" + FLAGS_output_tag;
  wtheta.Write(wtheta_file_name);
  wtheta.WriteCovariance(wcovar_file_name);
  //std::cout << "Writing galaxy cross-correlation to " <<
  //  wtheta_file_name << "\n";

  //std::ofstream output_file(wtheta_file_name.c_str());
  //for (Stomp::ThetaIterator iter=wtheta.Begin();iter!=wtheta.End();++iter) {
  //  output_file << std::setprecision(6) << iter->Theta() << " " <<
  //    iter->MeanWtheta()  << " " << iter->MeanWthetaError() << "\n";
  //}
  //output_file.close();

  return 0;
}
