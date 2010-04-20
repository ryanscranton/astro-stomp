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
DEFINE_string(input_file, "",
              "CSV list of ASCII file containing input circle parameters");
DEFINE_int32(max_resolution, 2048,
	     "Maximum resolution to use when pixelizing circles.");
DEFINE_int32(start_index, 0,
	     "Starting circle index.");
DEFINE_int32(finish_index, -1,
	     "Last circle index (default to all)");
DEFINE_double(weight, 1.0,
	      "Default weight for circle maps.");
DEFINE_bool(verbose, false,
	    "Print pixelization information for each circle individually.");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --output_file=<StompMap ASCII file>";
  usage += " --input_file=<ASCII circle parameter file>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_output_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_input_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  // First, we set up an empty Map to hold the pixelized circle Maps.
  Stomp::Map* stomp_map = new Stomp::Map();


  // Cast the maximum resolution to the proper variable type
  uint32_t max_resolution = static_cast<uint32_t>(FLAGS_max_resolution);

  // Now we extract the circle parameters, pixelize them and add their maps
  // to the global map.
  std::cout << "Parsing " << FLAGS_input_file << "...\n";
  std::ifstream circle_file(FLAGS_input_file.c_str());
  double ra, dec, radius;
  int32_t idx;

  if (!circle_file.is_open()) {
    std::cout << FLAGS_input_file << " does not exist!  Exiting.\n";
    exit(1);
  }

  int32_t n_circle = 0, n_kept = 0;
  double raw_area = 0.0, pixelized_raw_area = 0.0;
  uint32_t check = 1000;
  while (!circle_file.eof()) {
    circle_file >> idx >> ra >> dec >> radius;

    if (!circle_file.eof()) {
      if (FLAGS_start_index <= idx &&
	  (FLAGS_finish_index == -1 || FLAGS_finish_index > idx)) {
	Stomp::AngularCoordinate ang(ra, dec,
				     Stomp::AngularCoordinate::Equatorial);

	// Stars that are too close to the lambda poles won't be properly
	// pixelized, so we drop them.  They will be outside of our footprint
	// anyway, so there's no problem.
	if (Stomp::DoubleLE(ang.Lambda()+radius, 82.0) &&
	    Stomp::DoubleGE(ang.Lambda()-radius, -82.0)) {
	  Stomp::CircleBound* circle_bound =
	    new Stomp::CircleBound(ang, radius);
	  Stomp::Map* circle_map = new Stomp::Map(*circle_bound, 1.0,
					      max_resolution, FLAGS_verbose);
	  if (circle_map->Size() > 0) {
	    n_kept++;

	    pixelized_raw_area += circle_map->Area();
	    raw_area += circle_bound->Area();

	    stomp_map->IngestMap(*circle_map, true);

	    // If we're being verbose, output the results of pixelizing this
	    // circle.
	    if (FLAGS_verbose)
	      std::cout << "\t" << circle_map->Size() << " pixels, " <<
		circle_map->Area() << " sq. degrees. (" <<
		circle_bound->Area() << ")\n";
	  }
	  delete circle_map;
	  delete circle_bound;
	}
	n_circle++;
      }
      if (n_circle > 0 &&
	  (idx < FLAGS_finish_index || FLAGS_finish_index == -1) &&
	  n_circle % check == 0 && !FLAGS_verbose) {
	std::string status_file_name = FLAGS_output_file + "_status";
	std::ofstream status_file(status_file_name.c_str());
	status_file << idx << ": " << n_kept << "/" << n_circle << "/" <<
	  FLAGS_finish_index - FLAGS_start_index <<
	  " circles pixelized; " << pixelized_raw_area <<
	  " sq. degrees (" << raw_area << ")\n";
	status_file.close();
      }
    }
  }
  circle_file.close();

  if (stomp_map->Size() > 0) {
    std::cout << "Final map: " << stomp_map->Area() << " sq. degrees. (" <<
      pixelized_raw_area << " sq. degrees raw).\n" <<
      "Writing out to " << FLAGS_output_file << "...\n";

    stomp_map->Write(FLAGS_output_file);

    std::cout << "Done.\n";
  } else {
    std::cout << "No pixels in map.  Exiting.\n";
  }

  return 0;
}
