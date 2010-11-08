#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <stomp/stomp.h>
#include <gflags/gflags.h>


// Define our command-line flags.
DEFINE_string(map_file, "",
              "Name of the ASCII file containing the StompMap geometry");
DEFINE_string(sysmap_file, "",
              "Name of the ASCII file containing uniform resolution sysmap");
DEFINE_string(output_tag, "test",
              "Tag for output file: stripe_OUTPUT_TAG.hmap_combined");
DEFINE_double(seeing_min, 0.5, "Minimum allowable seeing value");
DEFINE_double(seeing_max, 2.0, "Maximum allowable seeing value");
DEFINE_double(extinction_min, 0.0, "Minimum allowable extinction value");
DEFINE_double(extinction_max, 1.0, "Maximum allowable extinction value");
DEFINE_double(sky_min, 0.0, "Minimum allowable sky brightness value");
DEFINE_double(sky_max, 24.0, "Maximum allowable sky brightness value");
DEFINE_bool(single_index, false, "Use older single-index file format.");
DEFINE_bool(no_weight, false, "Input map file is missing weight column.");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --map_file=<StompMap ASCII>";
  usage += " --sysmap_file=<SysMap ASCII>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_map_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_sysmap_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }


  // First, we read our STOMP map into a map object.  There are a couple
  // permutations based on the various map formats that are out there: with
  // or without a weight column or in the single index or double index format.
  std::cout << "Reading in initial map...\n";
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
  std::cout << "Read map from " << FLAGS_map_file << "; initial area: " <<
    stomp_map->Area() << " sq. deg.\n";


  // Now we read in our sysmap.  We only want to keep those pixels that
  // meet our selection criteria, so we filter on those as we read through
  // the file.
  Stomp::PixelVector sysmap_raw;
  std::ifstream sysmap_file(FLAGS_sysmap_file.c_str());
  double unmasked, seeing, extinction, sky;
  unsigned long hpixnum, superpixnum;
  unsigned int resolution, n_objects;

  unsigned long n_pixel = 0, n_keep = 0;
  while (!sysmap_file.eof()) {
    sysmap_file >> hpixnum >> superpixnum >> resolution >> unmasked >>
      seeing >> extinction >> sky >> n_objects;

    if (!sysmap_file.eof()) {
      n_pixel++;
      if ((seeing >= FLAGS_seeing_min) &&
	  (seeing <= FLAGS_seeing_max) &&
	  (extinction >= FLAGS_extinction_min) &&
	  (extinction <= FLAGS_extinction_max) &&
	  (sky >= FLAGS_sky_min) &&
	  (sky <= FLAGS_sky_max) &&
	  (n_objects > 0)) {
	Stomp::Pixel pix(resolution, hpixnum, superpixnum, 1.0);
	sysmap_raw.push_back(pix);
	n_keep++;
      }
    }
  }
  sysmap_file.close();
  std::cout << "Kept " << n_keep << "/" << n_pixel <<
    " systematics pixels; removing " <<
    sysmap_raw[0].Area()*(n_pixel - n_keep) << " sq. degrees...\n";


  // Convert that raw list of pixels into a StompMap.
  Stomp::Map* sysmap = new Stomp::Map(sysmap_raw);
  sysmap_raw.clear();


  // Now, we find the intersection between our allowed systematics map and
  // the basic map.
  std::cout << "Finding allowed area...\n";
  stomp_map->IntersectMap(*sysmap);


  // Finally write out the results...
  std::string output_file_name =
    "stripe_" + FLAGS_output_tag + ".hmap_combined";
  std::cout << "Writing map to " << output_file_name << ".  Final area: " <<
    stomp_map->Area() << " sq. deg.\n";

  stomp_map->Write(output_file_name);

  return 0;
}
