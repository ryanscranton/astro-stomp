#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <stomp.h>
#include <gflags/gflags.h>

// Define our command-line flags.
DEFINE_string(input_file, "",
              "Name of input ASCII file");
DEFINE_string(output_file, "",
              "Name of output ASCII file");
DEFINE_bool(survey_to_radec, false, "Convert from Lambda,Eta to Ra,Dec");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --input_file=<Input ASCII File>";
  usage += " --output_file=<Output ASCII File>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_input_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }
  if (FLAGS_output_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  Stomp::AngularCoordinate::Sphere coord_sphere = Stomp::AngularCoordinate::Equatorial;
  if (FLAGS_survey_to_radec) coord_sphere = Stomp::AngularCoordinate::Survey;

  std::ifstream coord_file(FLAGS_input_file.c_str());
  std::ofstream output_file(FLAGS_output_file.c_str());
  double theta, phi;

  while (!coord_file.eof()) {
    coord_file >> theta >> phi; 
    Stomp::AngularCoordinate tmp_ang(theta, phi,coord_sphere);
    if (!FLAGS_survey_to_radec) {
      output_file << tmp_ang.Lambda() << " " << tmp_ang.Eta() << "\n";
    }
    else {
      output_file << tmp_ang.RA() << " " << tmp_ang.DEC() << "\n";
    }
  }
  coord_file.close();
  output_file.close();

  return 0;
}
