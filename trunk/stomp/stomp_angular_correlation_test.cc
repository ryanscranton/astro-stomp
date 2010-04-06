#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>
#include "stomp_core.h"
#include "stomp_angular_bin.h"
#include "stomp_angular_correlation.h"
#include "stomp_map.h"
#include "stomp_scalar_map.h"

void AngularBinningTests() {
  // Now we break out the angular bin code.  This class lets you define either
  // a linear or logarithmic binning, depending on how you instantiate it.
  // We'll use the latter here.
  std::cout << "\n";
  std::cout << "****************************************\n";
  std::cout << "*** AngularCorrelation Binning Tests ***\n";
  std::cout << "****************************************\n";
  double theta = 3.0;
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);
  Stomp::ScalarMap* scalar_map =
    new Stomp::ScalarMap(*stomp_map, 128, Stomp::ScalarMap::DensityField);

  double theta_min = 0.01;
  double theta_max = 10.0;
  std::cout << "Begin by setting up a log-space binning from " <<
    theta_min << " to " << theta_max << "...\n";
  Stomp::AngularCorrelation *wtheta =
    new Stomp::AngularCorrelation(theta_min,theta_max,6.0,false);
  std::cout << "\t" << wtheta->NBins() << " angular bins:\n";
  for (Stomp::ThetaIterator iter=wtheta->Begin();iter!=wtheta->End();++iter)
    std::cout << "\t\t" << iter->ThetaMin() << " - " << iter->ThetaMax() <<
      " (" << iter->Sin2ThetaMin() << " - " << iter->Sin2ThetaMax() << ")\n";

  // Of course, just as important as the binning is knowing what resolution
  // map we need to use to calculate the auto-correlation on those scales.
  // This method figures that out.
  std::cout << "\tSetting up resolution values for each angular bin...\n";
  wtheta->AssignBinResolutions();
  for (Stomp::ThetaIterator iter=wtheta->Begin();iter!=wtheta->End();++iter)
    std::cout << "\t\t" << iter->Theta() << ": " << iter->Resolution() << "\n";

  // Now we test our methods for figuring out the angular bin extent we'd use
  // for a given map resolution.
  std::cout << "\tChecking to see which angular bins to use" <<
    " for our " << scalar_map->Resolution() << " density map...\n";
  std::cout << "\t\tAngular Range: " <<
    wtheta->ThetaMin(scalar_map->Resolution()) << " - " <<
    wtheta->ThetaMax(scalar_map->Resolution()) << "\n";

  std::cout << "\t\tIterator Test: " <<
    wtheta->Begin(scalar_map->Resolution())->ThetaMin() << " - " <<
    wtheta->End(scalar_map->Resolution())->ThetaMin() << "\n";

  // We've deviated a bit here from the STL to do the angular bin searching.
  // Instead of using the standard algorithms for figuring out which bin a
  // given angular separation (here given in terms of sin^2(theta) for reasons
  // which are clear if you look at the AutoCorrelate code), we use our own
  // binary search.  Given that, we need to make sure that it's doing the right
  // things.
  std::cout << "\tAngular Bin Search:\n";
  theta = 0.05;
  double sintheta =
    sin(theta*Stomp::DegToRad)*sin(theta*Stomp::DegToRad);
  Stomp::ThetaIterator theta_iter = wtheta->Find(wtheta->Begin(),wtheta->End(),
                                                 sintheta);
  std::cout << std::setprecision(6) <<
    "\t\t" << theta << ": " << theta_iter->ThetaMin() <<
    " - " << theta_iter->ThetaMax() << "\n";

  theta = 0.5;
  sintheta =
    sin(theta*Stomp::DegToRad)*sin(theta*Stomp::DegToRad);
  theta_iter = wtheta->Find(wtheta->Begin(),wtheta->End(),sintheta);
  std::cout << "\t\t" << theta << ": " << theta_iter->ThetaMin() <<
    " - " << theta_iter->ThetaMax() << "\n";

  theta = 5.0;
  sintheta =
    sin(theta*Stomp::DegToRad)*sin(theta*Stomp::DegToRad);
  theta_iter = wtheta->Find(wtheta->Begin(),wtheta->End(),sintheta);
  std::cout << std::setprecision(6) <<
    "\t\t" << theta << ": " << theta_iter->ThetaMin() <<
    " - " << theta_iter->ThetaMax() << "\n";

  // Now, restrict our search to just those bins that we'd use with our
  // density map and make sure that it's doing the right things.
  std::cout << "\tAngular Bin Search within the density map matching bins:\n";

  theta = 0.02;
  sintheta =
    sin(theta*Stomp::DegToRad)*sin(theta*Stomp::DegToRad);
  theta_iter = wtheta->Find(wtheta->Begin(scalar_map->Resolution()),
                            wtheta->End(scalar_map->Resolution()),sintheta);
  if (theta_iter == wtheta->End(scalar_map->Resolution())) {
    std::cout << "\t\tGood.\n" <<
      "\t\t\tTried to search for an angular bin outside the range" <<
      "\n\t\t\tand got back the end iterator.\n";
  } else {
    std::cout << "\t\tBad: " << theta << ": " << theta_iter->ThetaMin() <<
      " - " << theta_iter->ThetaMax() << "\n";
  }

  theta = 0.15;
  sintheta =
    sin(theta*Stomp::DegToRad)*sin(theta*Stomp::DegToRad);
  theta_iter = wtheta->Find(wtheta->Begin(scalar_map->Resolution()),
                            wtheta->End(scalar_map->Resolution()),sintheta);
  std::cout << "\t\t" << theta << ": " << theta_iter->ThetaMin() <<
    " - " << theta_iter->ThetaMax() << "\n";

  theta = 0.25;
  sintheta =
    sin(theta*Stomp::DegToRad)*sin(theta*Stomp::DegToRad);
  theta_iter = wtheta->Find(wtheta->Begin(scalar_map->Resolution()),
                            wtheta->End(scalar_map->Resolution()),sintheta);
  std::cout << "\t\t" << theta << ": " << theta_iter->ThetaMin() <<
    " - " << theta_iter->ThetaMax() << "\n";

  theta = 2.0;
  sintheta =
    sin(theta*Stomp::DegToRad)*sin(theta*Stomp::DegToRad);
  theta_iter = wtheta->Find(wtheta->Begin(scalar_map->Resolution()),
                            wtheta->End(scalar_map->Resolution()),sintheta);
  if (theta_iter == wtheta->End(scalar_map->Resolution())) {
    std::cout << "\t\tGood.\n" <<
      "\t\t\tTried to search for an angular bin outside the range" <<
      "\n\t\t\tand got back the end iterator.\n";
  } else {
    std::cout << "\t\tBad: " << theta << ": " << theta_iter->ThetaMin() <<
      " - " << theta_iter->ThetaMax() << "\n";
  }
}

// Define our command line flags
DEFINE_bool(all_angular_correlation_tests, false, "Run all class unit tests.");
DEFINE_bool(angular_binning_tests, false,
            "Run AngularCorrelation binning tests");

void AngularCorrelationUnitTests(bool run_all_tests) {
  void AngularBinningTests();

  if (run_all_tests) FLAGS_all_angular_correlation_tests = true;

  // Check the routines related to the AngularBin and AngularCorrelation
  // classes.
  if (FLAGS_all_angular_correlation_tests || FLAGS_angular_binning_tests)
    AngularBinningTests();
}
