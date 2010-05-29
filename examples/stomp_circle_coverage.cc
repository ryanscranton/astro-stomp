// Given a set of points and angular radii, we determine the fraction of the
// corresponding circular area that is contained in an input map.  We also
// perform corresponding checks on quadrants of the circle.
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <stomp.h>
#include <gflags/gflags.h>


// Define our command-line flags.
DEFINE_string(input_file, "",
              "Name of the ASCII file containing positions and radii");
DEFINE_string(map_file, "",
              "Name of the ASCII file containing the input map.");
DEFINE_string(output_file, "",
              "Name of the ASCII file to write results");
DEFINE_int32(max_resolution, -1,
             "Maximum resolution for pixelizing circles.  Default to map maximum resolution");
DEFINE_bool(single_index, false, "Use older single-index file input format.");
DEFINE_bool(no_weight, false, "Input map file is missing weight column.");
DEFINE_bool(radec, false, "Input positions are in equatorial coordinates.");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --input_file=<Input file with positions and radii>";
  usage += " --map_file=<StompMap ASCII>";
  usage += " --output_file=<Output file with quadrant data>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_input_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

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

  // First, we read our STOMP map into a map object.  There are a couple
  // permutations based on the various map formats that are out there: with
  // or without a weight column or in the single index or double index format.
  std::cout << "Reading in initial map...\n";
  Stomp::Map* stomp_map = new Stomp::Map(FLAGS_map_file, !FLAGS_single_index,
					 !FLAGS_no_weight);
  std::cout << "Read map from " << FLAGS_map_file << "; initial area: " <<
    stomp_map->Area() << " sq. deg.\n";

  // Now we scan through our input file and assemble our bitmask for each point.
  std::cout << "Parsing " << FLAGS_input_file << "...\n";
  std::ifstream input_file(FLAGS_input_file.c_str());
  std::ofstream output_file(FLAGS_output_file.c_str());

  Stomp::AngularCoordinate::Sphere coord = Stomp::AngularCoordinate::Survey;
  if (FLAGS_radec) coord = Stomp::AngularCoordinate::Equatorial;

  uint32_t max_resolution = stomp_map->MaxResolution();
  if (FLAGS_max_resolution != -1)
    max_resolution = static_cast<uint32_t>(FLAGS_max_resolution);
  std::cout << "Pixelizing at " << static_cast<int>(max_resolution) << "...\n";

  if (!input_file.is_open()) {
    std::cout << FLAGS_input_file << " does not exist!  Exiting.\n";
    exit(1);
  }

  // Define our bit mask values
  unsigned long inside_map = 1 << 0;
  unsigned long first_quadrant = 1 << 1;  // 0 <= position_angle < 90
  unsigned long second_quadrant = 1 << 2;  // 90 <= position_angle < 180
  unsigned long third_quadrant = 1 << 3;  // 180 <= position_angle < 270
  unsigned long fourth_quadrant = 1 << 4;  // 270 <= position_angle < 360

  int n_circle = 0;
  int check = 100;
  while (!input_file.eof()) {
    double theta, phi, radius;
    input_file >> theta >> phi >> radius;

    if (!input_file.eof()) {
      // Generate an AngularCoordinate based on our input position.
      Stomp::AngularCoordinate ang(theta, phi, coord);

      uint32_t idx = 0;

      // Figure out the first bit of our bitmask by checking that the position
      // is inside our map.
      if (stomp_map->Contains(ang)) idx += inside_map;

      // Now we assemble a Map for each of the quadrants and see if they are
      // inside our reference Map.
      {
	Stomp::WedgeBound wedge_bound(ang, radius, 0.0, 90.0, coord);
	Stomp::Map* wedge_map =
	  new Stomp::Map(wedge_bound, 1.0, max_resolution);

	if (stomp_map->Contains(*wedge_map)) idx += first_quadrant;

	delete wedge_map;
      }

      {
	Stomp::WedgeBound wedge_bound(ang, radius, 90.0, 180.0, coord);
	Stomp::Map* wedge_map =
	  new Stomp::Map(wedge_bound, 1.0, max_resolution);

	if (stomp_map->Contains(*wedge_map)) idx += second_quadrant;

	delete wedge_map;
      }

      {
	Stomp::WedgeBound wedge_bound(ang, radius, 180.0, 270.0, coord);
	Stomp::Map* wedge_map =
	  new Stomp::Map(wedge_bound, 1.0, max_resolution);

	if (stomp_map->Contains(*wedge_map)) idx += third_quadrant;

	delete wedge_map;
      }

      {
	Stomp::WedgeBound wedge_bound(ang, radius, 270.0, 360.0, coord);
	Stomp::Map* wedge_map =
	  new Stomp::Map(wedge_bound, 1.0, max_resolution);

	if (stomp_map->Contains(*wedge_map)) idx += fourth_quadrant;

	delete wedge_map;
      }

      n_circle++;
      if (n_circle > 0 && n_circle % check == 0)
	std::cout << "Processed " << n_circle << " objects...\n";

      output_file << theta << " " << phi << " " << radius << " " << idx << "\n";
    }
  }
  input_file.close();
  output_file.close();




  return 0;
}
