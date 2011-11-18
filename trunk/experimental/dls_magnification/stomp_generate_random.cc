#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <stomp.h>
#include <gflags/gflags.h>

// Define our command-line flags.
DEFINE_string(output_file, "",
              "Name of the ASCII file to store combined StompMap.");
DEFINE_string(input_map, "",
              "ASCII file containing the StompMap.");
DEFINE_string(unmasked_map, "",
	      "Unmasked fraction file to weight random points.");
DEFINE_string(region_file, "",
	      "ASCII file containing region definitions");
DEFINE_int32(n_randoms, 0, "Number of Random Points to Generate");
DEFINE_bool(galaxy_radec, true, "output Coordinates as RA, Dec");
DEFINE_bool(use_weights, false, "Generate Weighted random points");
DEFINE_bool(return_local_weight, false,
	    "Generate Uniform Random Points with weight from input_map");
DEFINE_bool(single_index, false, "Use older single-index file format.");
DEFINE_bool(no_weight, false, "Input file is missing weight column.");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --output_file=<Random Galaxy ASCII catalog>";
  usage += " --input_map=<StompMap ASCII file>";
  usage += " --n_randoms=<Number of Random Points to Generate>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_output_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_input_map.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }
  if (!(FLAGS_n_randoms>0)) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  // First we extract the input map names into a vector of strings.  The
  // input vector of file names should be comma-separated.
  Stomp::Map* stomp_map = new Stomp::Map(FLAGS_input_map.c_str());
  Stomp::Map* unmasked_map;
  if (!(FLAGS_unmasked_map.empty())) {
    unmasked_map = new Stomp::Map(FLAGS_unmasked_map.c_str());
  }

  std::cout << "Creating " << FLAGS_n_randoms << " random points \n"
            << "\ton Map: " << FLAGS_input_map << ", " 
	    << stomp_map->Area() << "sq. deg.\n" 
            << "\toutputing as " << FLAGS_output_file << "\n"; 

  Stomp::RegionBoundVector region_vector;
  if (!(FLAGS_region_file.empty())) {
    int type;
    int n_regions = 0;
    double region_area =0;
    double radius;
    double ra[4], dec[4];
    unsigned long idx;
    double weight;

    std::ifstream region_file(FLAGS_region_file.c_str());
    std::cout << "Reading regions from " << FLAGS_region_file; 
    while (!region_file.eof()) {
      region_file >> idx >> type >> radius >> weight 
		  >> ra[0] >> dec[0] >> ra[1] >> dec[1] 
		  >> ra[2] >> dec[2] >> ra[3] >> dec[3];
    
      if (!region_file.eof()) {
	
	Stomp::AngularVector ang_vector;
	for (int i=0;i<4;i++)
	  ang_vector.push_back(Stomp::AngularCoordinate(ra[i],dec[i],
			            Stomp::AngularCoordinate::Equatorial));
      
	Stomp::PolygonBound* poly_bound = new Stomp::PolygonBound(ang_vector);
	Stomp::RegionBound region_bound(poly_bound);
	region_area += region_bound.BoundArea();
	
	region_vector.push_back(region_bound);
	n_regions++;
      }
    }
    std::cout << "; " << n_regions << " Loaded; Total Area " << region_area 
	      << "; Average Area " << region_area/(1.0*n_regions) << "\n";
  }

  if (!FLAGS_return_local_weight) {
    Stomp::AngularVector random;
    stomp_map->GenerateRandomPoints(random, FLAGS_n_randoms, FLAGS_use_weights);
    
    std::ofstream output_file(FLAGS_output_file.c_str());
    for (Stomp::AngularIterator iter=random.begin();iter!=random.end();++iter) {
      output_file << iter->RA() << " " 
		  << iter->DEC() << " "
		  << 1.0 << " "
		  << 24.0 << "\n";
    }
    output_file.close();
  }
  else {

    double weight;
    std::ofstream output_file(FLAGS_output_file.c_str());
    long int n_randoms = 0;
    while (n_randoms < FLAGS_n_randoms) {
      Stomp::WeightedAngularCoordinate tmp_ang(0.0, 0.0, 0.0);
      stomp_map->GenerateSingleRandomPoint(tmp_ang, true, FLAGS_use_weights);

      if (!(FLAGS_unmasked_map.empty())) {
	weight = tmp_ang.Weight()*(unmasked_map->FindLocationWeight(tmp_ang));
      }
      else {
	weight = tmp_ang.Weight();
      }
      output_file << tmp_ang.RA() << " " 
		  << tmp_ang.DEC() << " ";
      if (!(FLAGS_region_file.empty())) {
	int region = 0;
	int hold_region = -1;
        double score = 1;
	double hold_score = 1;
	Stomp::Pixel tmp_pix = Stomp::Pixel(tmp_ang, 2048) ;
	for (uint32_t bound_iter=0;bound_iter<region_vector.size();
	     bound_iter++) {
	  score = region_vector[bound_iter].ScorePixel(tmp_pix);
	  if (score < hold_score) {
	    hold_score = score;
	    hold_region = region;
	  }
	  region++;
	}
	output_file << hold_region << " ";
      }
      output_file << weight << " "
		  << 24.0 << "\n";
      n_randoms++;
    }
    output_file.close();
  }

  return 0;
}
