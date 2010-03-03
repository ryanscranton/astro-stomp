#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>
#include "stomp_core.h"
#include "stomp_core_test.h"

void ConstantsTests() {
  // Check to make sure that our constants are all in place.
  std::cout << "\n";
  std::cout << "***********************\n";
  std::cout << "*** Constants Tests ***\n";
  std::cout << "***********************\n";
  std::cout << "\tpi: " << Stomp::Pi <<
    "\t\t\tdeg2Rad: " << Stomp::DegToRad << "\n";
  std::cout << "\trad2Deg: " << Stomp::RadToDeg <<
    "\t\tstrad2Deg: " << Stomp::StradToDeg << "\n";
  std::cout << "\tnx0: " << Stomp::Nx0 <<
    "\t\t\t\tny0: " << Stomp::Ny0 << "\n";
  std::cout << "\tetaOffSet: " << Stomp::EtaOffSet <<
    "\t\tsurveyCenterRA: " << Stomp::SurveyCenterRA << "\n";
  std::cout << "\tsurveyCenterDEC: " << Stomp::SurveyCenterDEC <<
    "\t\tnode: " << Stomp::Node << "\n";
  std::cout << "\tetaPole: " << Stomp::EtaPole <<
    "\t\thpix_resolution: " << Stomp::HPixResolution << "\n";
  std::cout << "\tmax_resolution: " << Stomp::MaxPixelResolution <<
    "\t\thpix_area: " << Stomp::HPixArea << "\n";
  std::cout << "\tmax_pixnum: " << Stomp::MaxPixnum <<
    "\t\tmax_superpixnum: " << Stomp::MaxSuperpixnum << "\n";
}

void MSBTests() {
  // Check to make sure that our code for finding ln2(int) is working
  std::cout << "\n";
  std::cout << "**********************************\n";
  std::cout << "*** Most Significant Bit Tests ***\n";
  std::cout << "**********************************\n";
  std::cout << "\t " << 4 << ": " <<
    static_cast<int>(Stomp::MostSignificantBit(4)) << "\n";
  std::cout << "\t " << 8 << ": " <<
    static_cast<int>(Stomp::MostSignificantBit(8)) << "\n";
  std::cout << "\t " << 16 << ": " <<
    static_cast<int>(Stomp::MostSignificantBit(16)) << "\n";
  std::cout << "\t " << 32 << ": " <<
    static_cast<int>(Stomp::MostSignificantBit(32)) << "\n";
  std::cout << "\t " << 64 << ": " <<
    static_cast<int>(Stomp::MostSignificantBit(64)) << "\n";
  std::cout << "\t " << 128 << ": " <<
    static_cast<int>(Stomp::MostSignificantBit(128)) << "\n";
  std::cout << "\t " << 256 << ": " <<
    static_cast<int>(Stomp::MostSignificantBit(256)) << "\n";
  std::cout << "\t " << 512 << ": " <<
    static_cast<int>(Stomp::MostSignificantBit(512)) << "\n";
  std::cout << "\t " << 1024 << ": " <<
    static_cast<int>(Stomp::MostSignificantBit(1024)) << "\n";
  std::cout << "\t " << 2048 << ": " <<
    static_cast<int>(Stomp::MostSignificantBit(2048)) << "\n";
  std::cout << "\t " << 4096 << ": " <<
    static_cast<int>(Stomp::MostSignificantBit(4096)) << "\n";
  std::cout << "\t " << 8192 << ": " <<
    static_cast<int>(Stomp::MostSignificantBit(8192)) << "\n";
  std::cout << "\t " << 16384 << ": " <<
    static_cast<int>(Stomp::MostSignificantBit(16384)) << "\n";
  std::cout << "\t " << 32768 << ": " <<
    static_cast<int>(Stomp::MostSignificantBit(32768)) << "\n";
}

int main(int argc, char **argv) {
  void ConstantsTests();
  void MSBTests();

  std::string usage = "Usage: ";
  usage += argv[0];
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  // First things first.  Let's check to make sure that our constants are
  // all in place.
  if (FLAGS_all_core_tests || FLAGS_core_constants_tests) ConstantsTests();

  // Now, we check our MostSignificantBit method to make sure it's working.
  if (FLAGS_all_core_tests || FLAGS_core_msb_tests) MSBTests();

  return 0;
}
