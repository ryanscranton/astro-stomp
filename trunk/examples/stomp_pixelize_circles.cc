#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include "stomp_util.h"
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


  // Now we extract the circle parameters, pixelize them and add their maps
  // to the global map.
  std::cout << "Parsing " << FLAGS_input_file << "...\n";
  std::ifstream circle_file(FLAGS_input_file.c_str());
  double ra, dec, radius;
  unsigned long idx;

  if (!circle_file.is_open()) {
    std::cout << FLAGS_input_file << " does not exist!  Exiting.\n";
    exit(1);
  }

  unsigned long n_circle = 0;
  unsigned long n_kept = 0;
  double raw_area = 0.0;
  double pixelized_raw_area = 0.0;
  unsigned long check = 1000;
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
	if (Stomp::Stomp::DoubleLE(ang.Lambda()+radius, 82.0) &&
	    Stomp::Stomp::DoubleGE(ang.Lambda()-radius, -82.0)) {
	  Stomp::CircleBound* circle_bound =
	    new Stomp::CircleBound(ang, radius, FLAGS_weight);

	  // Set the maximum resolution
	  circle_bound->SetMaxResolution(FLAGS_max_resolution);

	  // If we're being verbose, then we output the starting pixelization
	  // parameters.
	  if (FLAGS_verbose) {
	    int starting_resolution = circle_bound->FindStartingResolution();
	    circle_bound->FindXYBounds(starting_resolution);
	    std::cout << idx << ", (" << ang.Lambda() << "," << ang.Eta() <<
	      ", " << radius << "): " << starting_resolution << ", " <<
	      circle_bound->XMin() << " - " << circle_bound->XMax() << ", " <<
	      circle_bound->YMin() << " - " << circle_bound->YMax() << "\n";
	  }

	  // If the circle is too small to properly pixelize with the input
	  // resolution, then the pixelization may fail.  Provided that it
	  // doesn't, we can go forward and add the resulting Map to the
	  // aggregate.
	  if (circle_bound->Pixelize()) {
	    n_kept++;

	    pixelized_raw_area += circle_bound->PixelizedArea();
	    raw_area += circle_bound->Area();

	    Stomp::Map* circle_map = circle_bound->ExportMap();

	    stomp_map->IngestMap(*circle_map, true);

	    // If we're being verbose, output the results of pixelizing this
	    // circle.
	    if (FLAGS_verbose)
	      std::cout << "\t" << circle_bound->NPixel() <<
		" pixels, " << circle_bound->PixelizedArea() <<
		" sq. degrees. (" << circle_bound->Area() << ")\n";
	    delete circle_map;
	    delete circle_bound;
	  }
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
