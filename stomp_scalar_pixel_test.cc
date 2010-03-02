#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>
#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "stomp_pixel.h"
#include "stomp_scalar_pixel.h"

void StompScalarPixelTests() {
  // Before moving on to the map tests, we need to verify that the
  // Stomp::ScalarPixel, a derived class from Stomp::Pixel, is working
  // properly. Let's start with some initialization.
  std::cout << "\n";
  std::cout << "******************************\n";
  std::cout << "*** StompScalarPixel Tests ***\n";
  std::cout << "******************************\n";
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::ScalarPixel* tmp_scalar =
    new Stomp::ScalarPixel(ang, 256, 1.0, 1.0);
  Stomp::Pixel tmp_pix(ang, 256, 1.0);
  std::cout << "\tAngular Initialization: " <<
    tmp_scalar->Resolution() << " (" << tmp_pix.Resolution() << "), " <<
    tmp_scalar->PixelX() << " (" << tmp_pix.PixelX() << "), " <<
    tmp_scalar->PixelY() << " (" << tmp_pix.PixelY() << "), " <<
    tmp_scalar->Weight() << " (" << tmp_pix.Weight() << "), " <<
    tmp_scalar->Intensity() << "\n";

  delete tmp_scalar;

  tmp_scalar = new Stomp::ScalarPixel(tmp_pix.PixelX(), tmp_pix.PixelY(),
				      tmp_pix.Resolution(), 1.0, 1.0);
  std::cout << "\tIndex Initialization: " <<
    tmp_scalar->Resolution() << " (" << tmp_pix.Resolution() << "), " <<
    tmp_scalar->PixelX() << " (" << tmp_pix.PixelX() << "), " <<
    tmp_scalar->PixelY() << " (" << tmp_pix.PixelY() << "), " <<
    tmp_scalar->Weight() << " (" << tmp_pix.Weight() << "), " <<
    tmp_scalar->Intensity() << "\n";
  delete tmp_scalar;

  tmp_scalar = new Stomp::ScalarPixel(tmp_pix.Resolution(), tmp_pix.Pixnum(),
					1.0, 1.0);
  std::cout << "\tPixnum Initialization: " <<
    tmp_scalar->Resolution() << " (" << tmp_pix.Resolution() << "), " <<
    tmp_scalar->PixelX() << " (" << tmp_pix.PixelX() << "), " <<
    tmp_scalar->PixelY() << " (" << tmp_pix.PixelY() << "), " <<
    tmp_scalar->Weight() << " (" << tmp_pix.Weight() << "), " <<
    tmp_scalar->Intensity() << "\n";
  delete tmp_scalar;

  tmp_scalar = new Stomp::ScalarPixel(tmp_pix.Resolution(), tmp_pix.HPixnum(),
					tmp_pix.Superpixnum(), 1.0, 1.0);
  std::cout << "\tHPixnum Initialization: " <<
    tmp_scalar->Resolution() << " (" << tmp_pix.Resolution() << "), " <<
    tmp_scalar->PixelX() << " (" << tmp_pix.PixelX() << "), " <<
    tmp_scalar->PixelY() << " (" << tmp_pix.PixelY() << "), " <<
    tmp_scalar->Weight() << " (" << tmp_pix.Weight() << "), " <<
    tmp_scalar->Intensity() << "\n";
  delete tmp_scalar;

  tmp_scalar = new Stomp::ScalarPixel(ang, 32768, 1.0, 1.0);
  std::cout << "\tSpherical Coordinate Tests:\n";
  std::cout <<
    "\t\tLambda: " << tmp_scalar->Lambda() << " (" << ang.Lambda() <<
    "), Eta: " << tmp_scalar->Eta() << " (" << ang.Eta() << ")\n";
  std::cout << "\t\tX: " << tmp_scalar->UnitSphereX() <<
    " (" << ang.UnitSphereX() << ")\n";
  std::cout << "\t\tY: " << tmp_scalar->UnitSphereY() <<
    " (" << ang.UnitSphereY() << ")\n";
  std::cout << "\t\tZ: " << tmp_scalar->UnitSphereZ() <<
    " (" << ang.UnitSphereZ() << ")\n";
  std::cout <<
    "\t\tRA: " << tmp_scalar->RA() << " (" << ang.RA() <<
    "), DEC: " << tmp_scalar->DEC() << " (" << ang.DEC() << ")\n";
  std::cout <<
    "\t\tGalLat: " << tmp_scalar->GalLat() << " (" << ang.GalLat() <<
    "), GalLon: " << tmp_scalar->GalLon() << " (" << ang.GalLon() << ")\n";
}

DEFINE_bool(all_tests, false, "Run all unit tests.");
DEFINE_bool(scalar_pixel_tests, false, "Run ScalarPixel tests");

int main(int argc, char **argv) {
  void StompScalarPixelTests();

  std::string usage = "Usage: ";
  usage += argv[0];
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  // Check the basic inheritance of the Stomp::ScalarPixel class from
  // Stomp::Pixel.
  if (FLAGS_all_tests || FLAGS_scalar_pixel_tests) StompScalarPixelTests();

  return 0;
}
