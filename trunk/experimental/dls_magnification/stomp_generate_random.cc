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

  std::cout << "Creating " << FLAGS_n_randoms << " random points \n"
            << "\ton Map: " << FLAGS_input_map << ", " << stomp_map->Area() << "sq. deg.\n" 
            << "\toutputing as " << FLAGS_output_file << "\n"; 

  //if (FLAGS_single_index) {
  //  if (FLAGS_no_weight) {
  //    stomp_map = new Stomp::Map(FLAGS_input_map.c_str(), false, false);
  //  } else {
  //    stomp_map = new Stomp::Map(FLAGS_input_map.c_str(), false);
  //  }
  //} else {
  //  if (FLAGS_no_weight) {
  //    stomp_map = new Stomp::Map(FLAGS_input_map.c_str(), true, false);
  //  } else {
  //    stomp_map = new Stomp::Map(FLAGS_input_map.c_str());
  //  }
  //}

  if (!FLAGS_return_local_weight) {
    Stomp::AngularVector random;
    stomp_map->GenerateRandomPoints(random, FLAGS_n_randoms, FLAGS_use_weights);
    
    std::ofstream output_file(FLAGS_output_file.c_str());
    //output_file << "#fiat 1.0" << "\n"
    //<< "#ttype1 = raR" << "\n"
    //<< "#ttype2 = decR" << "\n"
    //<< "#ttype3 = weight" << "\n"
    //<< "##type4 = mag" << "\n";
    for (Stomp::AngularIterator iter=random.begin();iter!=random.end();++iter) {
      output_file << iter->RA() << " " 
		  << iter->DEC() << " "
		  << 1.0 << " "
		  << 24.0 << "\n";
    }
    output_file.close();
  }
  else {
    Stomp::WAngularVector random;
    while (random.size() < FLAGS_n_randoms) {
      Stomp::WeightedAngularCoordinate tmp_ang(0.0, 0.0, 0.0);
      stomp_map->GenerateSingleRandomPoint(tmp_ang, true, FLAGS_use_weights);
      random.push_back(tmp_ang);
    }

    std::ofstream output_file(FLAGS_output_file.c_str());
    //output_file << "#fiat 1.0" << "\n"
    //<< "#ttype1 = raR" << "\n"
    //<< "#ttype2 = decR" << "\n"
    //<< "#ttype3 = weight" << "\n"
    //<< "##type4 = mag" << "\n";
    for (Stomp::WAngularIterator iter=random.begin();
	 iter!=random.end();++iter) {
      output_file << iter->RA() << " " 
		  << iter->DEC() << " "
		  << iter->Weight() << " "
		  << 24.0 << "\n";
    }
    output_file.close();
  }

  return 0;
}
