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
DEFINE_string(output_file, "", 
              "Name of the ASCII file for the resultant StompMap");
DEFINE_string(stripes, "", "Comma-separated list of allowed stripes");
DEFINE_bool(single_index, false, "Use older single-index file format.");
DEFINE_bool(no_weight, false, "Input file is missing weight column.");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --map_file=<StompMap ASCII>";
  usage += " --output_file=<StompMap ASCII>";
  usage += " --stripes=<comma separated list of stripe numbers>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_map_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }
  if (FLAGS_output_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }
  if (FLAGS_stripes.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }


  // First, we read our STOMP map into a map object.  There are a couple
  // permutations based on the various map formats that are out there: with
  // or without a weight column or in the single index or double index format.
  std::cout << "Reading in initial map...\n";
  Stomp::Map* stomp_map = new Stomp::Map(FLAGS_map_file, !FLAGS_single_index,
					 !FLAGS_no_weight);
  std::cout << "Read map from " << FLAGS_map_file << "; initial area: " <<
    stomp_map->Area() << " sq. deg.\n";


  // Now we read in our list of allowed stripes.  This should be a
  // comma-separated list of stripe indices.  We use the static Stomp::Tokenize
  // method to convert the list of stripes into a vector of strings.
  std::vector<std::string> stripe_str;
  Stomp::Tokenize(FLAGS_stripes, stripe_str, ",");


  // Next, we iterate over the list of strings and convert them to stripe
  // integers.
  std::vector<int> stripes;
  std::cout << "Allowed stripes:\n";
  for (int i=0;i<stripe_str.size();++i) {
    int stripe_int = atoi(stripe_str[i].c_str());
    stripes.push_back(stripe_int);
    std::cout << "\t" << stripe_int << "\n";
  }


  // We're going to be doing some filtering on this, so we need to make sure
  // that the stripes are properly sorted.
  std::sort(stripes.begin(), stripes.end());


  // Now, we assemble vector of StompPixels from the allowed stripes that we'll
  // turn into a StompMap.
  Stomp::PixelVector pix;
  for (unsigned long k=0;k<Stomp::MaxSuperpixnum;k++) {
    Stomp::Pixel tmp_pix(Stomp::HPixResolution, k, 0.0);
    if (std::binary_search(stripes.begin(), stripes.end(), tmp_pix.Stripe()))
      pix.push_back(tmp_pix);
  }
  

  // Convert that raw list of pixels into a StompMap.
  Stomp::Map* stripe_map = new Stomp::Map(pix);
  pix.clear();


  // Now, we find the intersection between our allowed stripe map and
  // the basic map.
  std::cout << "Finding allowed area...\n";
  stomp_map->IntersectMap(*stripe_map);


  // Finally write out the results...
  std::cout << "Writing map to " << FLAGS_output_file << ".  Final area: " <<
    stomp_map->Area() << " sq. deg.\n";

  stomp_map->Write(FLAGS_output_file);

  return 0;
}
