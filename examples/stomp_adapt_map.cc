// This code uses the Map and ScalarMap objects to attempt to match an input
// map to data.  The idea is that the input map is a close, but not quite
// accurate description of the area covered by the data.  By building a
// ScalarMap with the input Map and then populating it with points from the
// data file, we can get a reasonable approximation of which parts of the input
// Map are within the area spanned by the data file.

#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <stomp.h>
#include <gflags/gflags.h>

// Define our command-line flags.
DEFINE_string(input_map, "",
              "Name of the ASCII file containing the Map geometry");
DEFINE_string(galaxy_file, "",
              "Name of the ASCII file containing the input galaxy catalog");
DEFINE_string(output_map, "",
              "Name of the ASCII file target for writing the resulting Map");
DEFINE_bool(single_index, false, "Use older single-index file format.");
DEFINE_bool(no_weight, false, "Input file is missing weight column.");
DEFINE_int32(resolution, 128,
             "Resolution to use for the ScalarMap");
DEFINE_bool(radec, false,
	    "Input file is in equatorial coordinates (default is survey)");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --input_map=<StompMap ASCII>";
  usage += " --galaxy_file=<Galaxy catalog ASCII>";
  usage += " --output_map=<StompMap ASCII>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_input_map.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_galaxy_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_output_map.empty()) {
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
      stomp_map = new Stomp::Map(FLAGS_input_map, false, false);
    } else {
      stomp_map = new Stomp::Map(FLAGS_input_map, false);
    }
  } else {
    if (FLAGS_no_weight) {
      stomp_map = new Stomp::Map(FLAGS_input_map, true, false);
    } else {
      stomp_map = new Stomp::Map(FLAGS_input_map);
    }
  }
  std::cout << "Read map from " << FLAGS_input_map << "; total area: " <<
    stomp_map->Area() << " sq. deg.\n";

  // Now we create our ScalarMap from the Map at the specified resolution.
  uint32_t resolution  = static_cast<uint32_t>(FLAGS_resolution);
  Stomp::ScalarMap* scalar_map =
    new Stomp::ScalarMap(*stomp_map, resolution,
			 Stomp::ScalarMap::DensityField);

  // Now we read in our galaxy data file.  The expected format is
  //  LAMBDA  ETA  WEIGHT  MAGNITUDE
  // where the WEIGHT column is the likelihood that the object is a galaxy
  // and MAGNITUDE is the apparent magnitude in a given filter.  We filter all
  // of the objects against the map, tossing out any objects that aren't in the
  // map.
  std::ifstream galaxy_file(FLAGS_galaxy_file.c_str());
  double theta, phi;
  uint32_t n_galaxy = 0;

  Stomp::AngularCoordinate::Sphere coord = Stomp::AngularCoordinate::Survey;
  if (FLAGS_radec) coord = Stomp::AngularCoordinate::Equatorial;

  while (!galaxy_file.eof()) {
    galaxy_file >> theta >> phi;
    Stomp::AngularCoordinate tmp_ang(theta, phi, coord);

    if (stomp_map->FindLocation(tmp_ang)) {
      scalar_map->AddToMap(tmp_ang);
    }

    n_galaxy++;
  }
  galaxy_file.close();

  std::cout << "Read " << n_galaxy << " galaxies from " << FLAGS_galaxy_file <<
    "; kept " << scalar_map->NPoints() << "\n";

  // Now we iterate over our ScalarMap contents and copy the pixels that contain
  // data to a PixelVector.
  std::cout << "Finding filled pixels...\n";
  Stomp::PixelVector common_pix;
  for (Stomp::ScalarIterator iter=scalar_map->Begin();
       iter!=scalar_map->End();++iter) {
    if (iter->NPoints() > 0) {
      common_pix.push_back(Stomp::Pixel(iter->PixelX(), iter->PixelY(),
					iter->Resolution(), 1.0));
    }
  }

  // We're now done with the ScalarMap, so we can free up that memory.
  delete scalar_map;

  // Now we make a Map out of the Pixels that contained data points
  Stomp::Map* common_map = new Stomp::Map(common_pix);

  // And find the intersection between this Map and our original one.  At this
  // point we're done and can write the results to the output Map file.
  std::cout << "Calculating intersection...\n";
  if (stomp_map->IntersectMap(*common_map)) {
    std::cout << "Writing results to " << FLAGS_output_map << "...\n";
    stomp_map->Write(FLAGS_output_map);
  } else {
    std::cout << "Failed to find any intersecting area between the two Maps.\n";
  }

  return 0;
}
