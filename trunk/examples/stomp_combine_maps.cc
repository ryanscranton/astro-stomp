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
DEFINE_string(input_files, "",
              "CSV list of ASCII files containing input StompMaps.");
DEFINE_string(intersect_files, "",
              "CSV list of ASCII files containing StompMaps to intersect.");
DEFINE_string(exclude_files, "",
              "CSV list of ASCII files containing StompMaps to exclude.");
DEFINE_string(coverage_file, "",
              "Name of the ASCII file to store output Map Coverage.");
DEFINE_int32(max_resolution, -1,
	     "Maximum resolution to use for output map.");
DEFINE_double(max_weight, 1.0,
	      "Maximum allowed weight for output map.");
DEFINE_double(min_weight, 0.0,
	      "Minimum allowed weight for output map.");
DEFINE_bool(average_weights, false, "Average softened pixel weights.");
DEFINE_bool(single_index, false, "Use older single-index file format.");
DEFINE_bool(no_weight, false, "Input file is missing weight column.");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --output_file=<StompMap ASCII file>";
  usage += " --input_files=<CSV list of StompMap ASCII files>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_output_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_input_files.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  // First we extract the input map names into a vector of strings.  The
  // input vector of file names should be comma-separated.
  std::vector<std::string> input_files;
  Stomp::Tokenize(FLAGS_input_files, input_files, ",");
  Stomp::Map* stomp_map;
  double raw_area = 0.0;

  if (input_files.size() > 1) {
    std::cout << "Combining " << input_files.size() << " files...\n";

    stomp_map = new Stomp::Map();

    for (std::vector<std::string>::iterator file_iter=input_files.begin();
	 file_iter!=input_files.end();++file_iter) {

      std::cout << "\tParsing " << file_iter->c_str() << "...\n";

      Stomp::Map* tmp_map;
      if (FLAGS_single_index) {
	if (FLAGS_no_weight) {
	  tmp_map = new Stomp::Map(*file_iter, false, false);
	} else {
	  tmp_map = new Stomp::Map(*file_iter, false);
	}
      } else {
	if (FLAGS_no_weight) {
	  tmp_map = new Stomp::Map(*file_iter, true, false);
	} else {
	  tmp_map = new Stomp::Map(*file_iter);
	}
      }

      std::cout << "\t\tAdding map with " << tmp_map->Area() <<
	" sq. degrees (" << tmp_map->Size() << ")...\n";

      std::cout << "check\n";
      raw_area += tmp_map->Area();
      std::cout << "check\n";
      stomp_map->IngestMap(*tmp_map, true);
      std::cout << "check\n";
      tmp_map->Clear();
      std::cout << "check\n";
      delete tmp_map;
      std::cout << "check\n";
    }
  } else {
    std::cout << "\tParsing " << input_files[0].c_str() << "...\n";

    if (FLAGS_single_index) {
      if (FLAGS_no_weight) {
	stomp_map = new Stomp::Map(input_files[0], false, false);
      } else {
	stomp_map = new Stomp::Map(input_files[0], false);
      }
    } else {
      if (FLAGS_no_weight) {
	stomp_map = new Stomp::Map(input_files[0], true, false);
      } else {
	stomp_map = new Stomp::Map(input_files[0]);
      }
    }

    std::cout << "\t\tAdding map with " << stomp_map->Area() <<
      " sq. degrees...\n";
    raw_area = stomp_map->Area();
  }

  std::vector<std::string> intersect_files;
  Stomp::Tokenize(FLAGS_intersect_files, intersect_files, ",");
  if (intersect_files.size() > 0) {
    std::cout << "Finding intersection against " <<
      intersect_files.size() << " files...\n";

    for (std::vector<std::string>::iterator file_iter=intersect_files.begin();
	 file_iter!=intersect_files.end();++file_iter) {

      std::cout << "\tParsing " << file_iter->c_str() << "...\n";

      Stomp::Map* tmp_map;
      if (FLAGS_single_index) {
	if (FLAGS_no_weight) {
	  tmp_map = new Stomp::Map(*file_iter, false, false);
	} else {
	  tmp_map = new Stomp::Map(*file_iter, false);
	}
      } else {
	if (FLAGS_no_weight) {
	  tmp_map = new Stomp::Map(*file_iter, true, false);
	} else {
	  tmp_map = new Stomp::Map(*file_iter);
	}
      }

      std::cout << "\t\tFinding intersection with map with " <<
	tmp_map->Area() << " sq. degrees...\n";

      stomp_map->IntersectMap(*tmp_map);
      tmp_map->Clear();
      delete tmp_map;
    }
    raw_area = stomp_map->Area();
    std::cout << "\tCurrent area: " << stomp_map->Area() << " sq. degrees.\n";
  }

  std::vector<std::string> exclude_files;
  Stomp::Tokenize(FLAGS_exclude_files, exclude_files, ",");
  if (exclude_files.size() > 0) {
    std::cout << "Excluding area from " <<
      exclude_files.size() << " files...\n";

    for (std::vector<std::string>::iterator file_iter=exclude_files.begin();
	 file_iter!=exclude_files.end();++file_iter) {

      std::cout << "\tParsing " << file_iter->c_str() << "...\n";

      Stomp::Map* tmp_map;
      if (FLAGS_single_index) {
	if (FLAGS_no_weight) {
	  tmp_map = new Stomp::Map(*file_iter, false, false);
	} else {
	  tmp_map = new Stomp::Map(*file_iter, false);
	}
      } else {
	if (FLAGS_no_weight) {
	  tmp_map = new Stomp::Map(*file_iter, true, false);
	} else {
	  tmp_map = new Stomp::Map(*file_iter);
	}
      }

      std::cout << "\t\tExcluding map with " << tmp_map->Area() <<
	" sq. degrees...\n";

      if (stomp_map->ExcludeMap(*tmp_map, false)) {
	std::cout << "\t\t\tArea remaining: " << stomp_map->Area() <<
	  " sq. degrees...\n";
      } else {
	std::cout << "\t\t\tWARNING: No area remaining in map!\n";
      }
      tmp_map->Clear();
      delete tmp_map;
    }
    raw_area = stomp_map->Area();
    std::cout << "\tCurrent area: " << stomp_map->Area() << " sq. degrees.\n";
  }

  if ((FLAGS_max_resolution > 0) && (FLAGS_max_resolution % 2 == 0)) {
    std::cout << "Softening pixels with resolution above " <<
      FLAGS_max_resolution << "\n\t";
    if (FLAGS_average_weights) {
      std::cout << "Averaging weights...\n";
    } else {
      std::cout << "Using unity weights...\n";
    }
    std::cout << "\tStarting area: " << stomp_map->Area() << " sq. degrees.\n";
    stomp_map->Soften(FLAGS_max_resolution, FLAGS_average_weights);
    std::cout << "\tEnding area: " << stomp_map->Area() << " sq. degrees.\n";
  }

  std::cout << "Filtering weight: " << FLAGS_min_weight << " <= weight <= " <<
    FLAGS_max_weight << "\n";
  std::cout << "\tStarting area: " << stomp_map->Area() << " sq. degrees.\n";
  stomp_map->SetMinimumWeight(FLAGS_min_weight);
  stomp_map->SetMaximumWeight(FLAGS_max_weight);
  std::cout << "\tEnding area: " << stomp_map->Area() << " sq. degrees.\n";
  std::cout << "\t\t" << stomp_map->MinWeight() << " <= weight <= " <<
    stomp_map->MaxWeight() << "\n";


  std::cout << "Final map: " << stomp_map->Area() << " sq. degrees. (" <<
    raw_area << " sq. degrees raw).\n" <<
    "Writing out to " << FLAGS_output_file << "...\n";

  stomp_map->Write(FLAGS_output_file);

  std::cout << "Done.\n";

  return 0;
}
