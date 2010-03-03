#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>
#include "stomp_util.h"
#include "stomp_util_test.h"

void CosmologyTests() {
  // Check to make sure that our Cosmology static class works properly.
  std::cout << "\n";
  std::cout << "***********************\n";
  std::cout << "*** Cosmology Tests ***\n";
  std::cout << "***********************\n";

  // First check the visible constant values.
  std::cout << "Cosmology: Omega_M = " << Stomp::Cosmology::OmegaM() <<
    ", Omega_L = " << Stomp::Cosmology::OmegaL() << ", Hubble Constant = " <<
    Stomp::Cosmology::HubbleConstant() << "\n";

  // Now we check some distances.
  std::cout << "Comoving distances:\n";
  double z = 0.0;
  for (uint16_t i=0;i<20;i++) {
    std::cout << "\tD(" << z << ") = " <<
      Stomp::Cosmology::ComovingDistance(z) << " Mpc/h\n";
    z += 0.2;
  }

  // Now, change the cosmology slightly and re-calculate.
  Stomp::Cosmology::SetOmegaM(0.3);
  std::cout << "\nComoving distances for Omega_M = " <<
    Stomp::Cosmology::OmegaM() << ":\n";
  z = 0.0;
  for (uint16_t i=0;i<20;i++) {
    std::cout << "\tD(" << z << ") = " <<
      Stomp::Cosmology::ComovingDistance(z) << " Mpc/h\n";
    z += 0.2;
  }
}

void UtilUnitTests(bool run_all_tests) {
  void StompCosmologyTests();

  if (run_all_tests) FLAGS_all_util_tests = true;

  // Now, we check our static Cosmology class to make sure it's functioning.
  if (FLAGS_all_util_tests || FLAGS_util_cosmology_tests) CosmologyTests();
}
