#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>
#include "stomp_core.h"
#include "stomp_util.h"
#include "stomp_angular_coordinate.h"
#include "stomp_pixel.h"
#include "stomp_pixel_test.h"

void PixelBasicTests() {
  // Ok, now we'll try moving around a bit in resolution space.
  std::cout << "\n";
  std::cout << "***************************\n";
  std::cout << "*** Pixel Basic scaling ***\n";
  std::cout << "***************************\n";
  std::cout << "Setting index at each step manually:\n";

  Stomp::AngularCoordinate ang(20.0, 20.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);

  for (uint16_t resolution=Stomp::MaxPixelResolution;
       resolution>=Stomp::HPixResolution;resolution/=2) {
    tmp_pix.SetResolution(resolution);
    tmp_pix.SetPixnumFromAng(ang);
    std::cout << "\tHPixnum, Superpixnum, Resolution: " <<
      tmp_pix.HPixnum() << ", " << tmp_pix.Superpixnum() <<
      ", " << tmp_pix.Resolution() << "\n";
  }
  std::cout << "Now scaling with the Superpix function:\n";
  tmp_pix.SetResolution(Stomp::MaxPixelResolution);
  tmp_pix.SetPixnumFromAng(ang);

  for (uint16_t resolution=Stomp::MaxPixelResolution;
       resolution>=Stomp::HPixResolution;resolution/=2) {
    tmp_pix.SetToSuperPix(resolution);
    std::cout << "\tHPixnum, Superpixnum, Resolution: " <<
      tmp_pix.HPixnum() << ", " << tmp_pix.Superpixnum() <<
      ", " << tmp_pix.Resolution() << "\n";
  }

  //  Ok, now we go the other way.
  std::cout << "\nSub-pixel index checking:\n";

  tmp_pix.SetResolution(Stomp::HPixResolution);
  tmp_pix.SetPixnumFromAng(ang);
  std::cout << "\tX, Y, Resolution: " <<
    tmp_pix.PixelX() << ", " << tmp_pix.PixelY() <<
    ", " << tmp_pix.Resolution() << "\n";

  for (uint16_t resolution=Stomp::HPixResolution, i=0;
       i<Stomp::ResolutionLevels;resolution*=2,i++) {
    uint32_t x_min, x_max, y_min, y_max, n_pixel;
    tmp_pix.SubPix(resolution,x_min,x_max,y_min,y_max);
    n_pixel = (x_max - x_min + 1)*(y_max - y_min + 1);
    std::cout << "\t" << resolution <<
      ", X: " << x_min << " - " << x_max <<
      ", Y: " << y_min << " - " << y_max <<
      ", " << n_pixel << " pixels\n";
  }
}

void PixelStripeTests() {
  // Now, we'll do some spherical geometry Tests.  First up, we check the
  // multi-resolution stripe function.
  std::cout << "\n";
  std::cout << "*************************\n";
  std::cout << "*** Pixel Stripe test ***\n";
  std::cout << "*************************\n";
  Stomp::AngularCoordinate ang(20.0, 20.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);
  tmp_pix.SetResolution(Stomp::MaxPixelResolution);
  tmp_pix.SetPixnumFromAng(ang);

  for (uint16_t resolution=Stomp::MaxPixelResolution;
       resolution>=Stomp::HPixResolution;resolution/=2) {
    std::cout << "\tResolution " << resolution << ": Stripe = " <<
      tmp_pix.Stripe(resolution) << "\n";
  }
}

void PixelXYTests() {
  // Now, the XY bounds for a given angular range around a given point.  There
  // are two things to note here.  First, finding all of the pixels within
  // a given angular radius means looking at many more pixels at high
  // latitude than at low latitude.  Second, if we use the flexible box version
  // of XYBounds (where the x bounds are calculated separately for each value
  // of y), then the number of pixels that we'd examine is much smaller at
  // high latitudes.  Hence, even though this version of XYBounds is slower,
  // you may save time (and memory) on the back end if you have to do
  // several operations on the pixels you get back.
  std::cout << "\n";
  std::cout << "***********************\n";
  std::cout << "*** X-Y bounds test ***\n";
  std::cout << "***********************\n";
  Stomp::AngularCoordinate ang(0.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);
  double theta = 10.0;  //  Look for all pixels within 10 degrees.

  std::cout << "Low latitude (Lambda = 0.0), fixed box\n";
  ang.SetSurveyCoordinates(0.0,0.0);
  tmp_pix.SetResolution(Stomp::MaxPixelResolution);
  tmp_pix.SetPixnumFromAng(ang);
  for (uint16_t resolution=Stomp::MaxPixelResolution;
       resolution>=Stomp::HPixResolution;resolution/=2) {
    tmp_pix.SetToSuperPix(resolution);

    uint32_t x_min, x_max, y_min, y_max, n_pixel;
    tmp_pix.XYBounds(theta,x_min,x_max,y_min,y_max,true);
    n_pixel = (x_max - x_min + 1)*(y_max - y_min + 1);
    std::cout << "\t" << resolution <<
      ", X: " << x_min << " - " << x_max <<
      ", Y: " << y_min << " - " << y_max <<
      ", " << n_pixel << " pixels\n";
  }

  std::cout << "High latitude (Lambda = 60.0), fixed box\n";

  ang.SetSurveyCoordinates(60.0,0.0);
  tmp_pix.SetResolution(Stomp::MaxPixelResolution);
  tmp_pix.SetPixnumFromAng(ang);
  for (uint16_t resolution=Stomp::MaxPixelResolution;
       resolution>=Stomp::HPixResolution;resolution/=2) {
    tmp_pix.SetToSuperPix(resolution);

    uint32_t x_min, x_max, y_min, y_max, n_pixel;
    tmp_pix.XYBounds(theta,x_min,x_max,y_min,y_max,true);
    n_pixel = (x_max - x_min + 1)*(y_max - y_min + 1);
    std::cout << "\t" << resolution <<
      ", X: " << x_min << " - " << x_max <<
      ", Y: " << y_min << " - " << y_max <<
      ", " << n_pixel << " pixels\n";
  }

  std::cout << "Low latitude (Lambda = 0.0), flexible box\n";

  ang.SetSurveyCoordinates(0.0,0.0);
  tmp_pix.SetResolution(Stomp::MaxPixelResolution);
  tmp_pix.SetPixnumFromAng(ang);
  for (uint16_t resolution=Stomp::MaxPixelResolution;
       resolution>=Stomp::HPixResolution;resolution/=2) {
    tmp_pix.SetToSuperPix(resolution);

    uint32_t y_min, y_max, n_pixel;
    std::vector<uint32_t> x_min, x_max;
    tmp_pix.XYBounds(theta,x_min,x_max,y_min,y_max,true);
    n_pixel = 0;
    for (uint32_t y=y_min,m=0;y<=y_max;y++,m++)
      n_pixel += x_max[m] - x_min[m] + 1;
    std::cout << "\t" << resolution <<
      ", Y: " << y_min << " - " << y_max <<
      ", " << n_pixel << " pixels\n";
  }

  std::cout << "High latitude (Lambda = 60.0), flexible box\n";

  ang.SetSurveyCoordinates(60.0,0.0);
  tmp_pix.SetResolution(Stomp::MaxPixelResolution);
  tmp_pix.SetPixnumFromAng(ang);
  for (uint16_t resolution=Stomp::MaxPixelResolution;
       resolution>=Stomp::HPixResolution;resolution/=2) {
    tmp_pix.SetToSuperPix(resolution);

    uint32_t y_min, y_max, n_pixel;
    std::vector<uint32_t> x_min, x_max;
    tmp_pix.XYBounds(theta,x_min,x_max,y_min,y_max,true);
    n_pixel = 0;
    for (uint32_t y=y_min,m=0;y<=y_max;y++,m++)
      n_pixel += x_max[m] - x_min[m] + 1;
    std::cout << "\t" << resolution <<
      ", Y: " << y_min << " - " << y_max <<
      ", " << n_pixel << " pixels\n";
  }
}

void PixelBoundTests() {
  // Quick checks for the various routines to return the pixel bounds in the
  // various coordinate systems.
  std::cout << "\n";
  std::cout << "****************************\n";
  std::cout << "*** Pixel Boundary Tests ***\n";
  std::cout << "****************************\n";

  uint16_t starting_resolution = 2048;
  Stomp::AngularCoordinate ang(0.0, 180.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, starting_resolution);

  // First we do things for Survey coordinates at the Eta discontinuity.  This
  // should be the simplest case since this is the natural coordinate system
  // for the pixels.
  std::cout << "\nSurvey Coordinates (Lambda = " << ang.Lambda() <<
    ", Eta = " << ang.Eta() << "):\n";
  for (uint16_t resolution=starting_resolution;
       resolution>=Stomp::HPixResolution;resolution /= 2) {
    tmp_pix.SetToSuperPix(resolution);
    std::cout << "\t" << tmp_pix.Resolution() << ": " << tmp_pix.EtaMin() <<
      " < Eta < " << tmp_pix.EtaMax() << " (" <<
      tmp_pix.EtaMaxContinuous() << " continuous)\n";
  }

  // Now we do the same for the Equatorial coordinates at RA = 0, DEC = 0.
  ang.SetEquatorialCoordinates(0.0, 0.0);
  tmp_pix.SetResolution(starting_resolution);
  tmp_pix.SetPixnumFromAng(ang);
  std::cout << "\nEquatorial Coordinates (RA = " << ang.RA() <<
    ", DEC = " << ang.DEC() << "):\n";
  for (uint16_t resolution=starting_resolution;
       resolution>=Stomp::HPixResolution;resolution /= 2) {
    tmp_pix.SetToSuperPix(resolution);
    std::cout << "\t" << tmp_pix.Resolution() << ": " << tmp_pix.RAMin() <<
      " < RA < " << tmp_pix.RAMax() << " (" <<
      tmp_pix.RAMaxContinuous() << " continuous)\n";
  }

  // Now we move to a higher DEC to check things there.
  ang.SetEquatorialCoordinates(0.0, 75.0);
  tmp_pix.SetResolution(starting_resolution);
  tmp_pix.SetPixnumFromAng(ang);
  std::cout << "\nEquatorial Coordinates (RA = " << ang.RA() <<
    ", DEC = " << ang.DEC() << "):\n";
  for (uint16_t resolution=starting_resolution;
       resolution>=Stomp::HPixResolution;resolution /= 2) {
    tmp_pix.SetToSuperPix(resolution);
    std::cout << "\t" << tmp_pix.Resolution() << ": " << tmp_pix.RAMin() <<
      " < RA < " << tmp_pix.RAMax() << " (" <<
      tmp_pix.RAMaxContinuous() << " continuous)\n";
  }

  // And a final check near the pole.
  ang.SetEquatorialCoordinates(0.0, 88.0);
  tmp_pix.SetResolution(starting_resolution);
  tmp_pix.SetPixnumFromAng(ang);
  std::cout << "\nEquatorial Coordinates (RA = " << ang.RA() <<
    ", DEC = " << ang.DEC() << "):\n";
  for (uint16_t resolution=starting_resolution;
       resolution>=Stomp::HPixResolution;resolution /= 2) {
    tmp_pix.SetToSuperPix(resolution);
    std::cout << "\t" << tmp_pix.Resolution() << ": " << tmp_pix.RAMin() <<
      " < RA < " << tmp_pix.RAMax() << " (" <<
      tmp_pix.RAMaxContinuous() << " continuous)\n";
  }

  // Now we do the same for the Galactic coordinates at Lon = 0, Lat = 0.
  ang.SetGalacticCoordinates(0.0, 0.0);
  tmp_pix.SetResolution(starting_resolution);
  tmp_pix.SetPixnumFromAng(ang);
  std::cout << "\nGalactic Coordinates (GalLon = " << ang.GalLon() <<
    ", GalLat = " << ang.GalLat() << "):\n";
  for (uint16_t resolution=starting_resolution;
       resolution>=Stomp::HPixResolution;resolution /= 2) {
    tmp_pix.SetToSuperPix(resolution);
    std::cout << "\t" << tmp_pix.Resolution() << ": " << tmp_pix.GalLonMin() <<
      " < GalLon < " << tmp_pix.GalLonMax() << " (" <<
      tmp_pix.GalLonMaxContinuous() << " continuous)\n";
  }

  // Now something at higher latitude.
  ang.SetGalacticCoordinates(0.0, 75.0);
  tmp_pix.SetResolution(starting_resolution);
  tmp_pix.SetPixnumFromAng(ang);
  std::cout << "\nGalactic Coordinates (GalLon = " << ang.GalLon() <<
    ", GalLat = " << ang.GalLat() << "):\n";
  for (uint16_t resolution=starting_resolution;
       resolution>=Stomp::HPixResolution;resolution /= 2) {
    tmp_pix.SetToSuperPix(resolution);
    std::cout << "\t" << tmp_pix.Resolution() << ": " << tmp_pix.GalLonMin() <<
      " < GalLon < " << tmp_pix.GalLonMax() << " (" <<
      tmp_pix.GalLonMaxContinuous() << " continuous)\n";
  }

  // And a final check near the pole.
  ang.SetGalacticCoordinates(0.0, 88.0);
  tmp_pix.SetResolution(starting_resolution);
  tmp_pix.SetPixnumFromAng(ang);
  std::cout << "\nGalactic Coordinates (GalLon = " << ang.GalLon() <<
    ", GalLat = " << ang.GalLat() << "):\n";
  for (uint16_t resolution=starting_resolution;
       resolution>=Stomp::HPixResolution;resolution /= 2) {
    tmp_pix.SetToSuperPix(resolution);
    std::cout << "\t" << tmp_pix.Resolution() << ": " << tmp_pix.GalLonMin() <<
      " < GalLon < " << tmp_pix.GalLonMax() << " (" <<
      tmp_pix.GalLonMaxContinuous() << " continuous)\n";
  }
}

void PixelWithinRadiusTests() {
  // Now we check our WithinRadius routine.  We'll stick with a 10 degree
  // radius, but cut the maximum resolution to check to 256 to save time.
  std::cout << "\n";
  std::cout << "******************************************************\n";
  std::cout << "*** Returning pixels within an angular radius test ***\n";
  std::cout << "******************************************************\n";

  Stomp::AngularCoordinate ang(0.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 128);
  double theta = 10.0;  //  Look for all pixels within 10 degrees.
  std::cout << "Low latitude (Lambda = 0.0), flexible box\n";
  for (uint16_t resolution=128;
       resolution>=Stomp::HPixResolution;resolution/=2) {
    tmp_pix.SetToSuperPix(resolution);

    uint32_t y_min, y_max, n_pixel;
    std::vector<uint32_t> x_min, x_max;
    tmp_pix.XYBounds(theta,x_min,x_max,y_min,y_max,true);
    n_pixel = 0;
    for (uint32_t y=y_min,m=0;y<=y_max;y++,m++)
      n_pixel += x_max[m] - x_min[m] + 1;

    Stomp::PixelVector annulus_pix;
    tmp_pix.WithinRadius(theta,annulus_pix);
    std::cout << "\t" << resolution <<
        ", Checked " << n_pixel << " pixels, kept " <<
        annulus_pix.size() << "\n";
  }

  std::cout << "High latitude (Lambda = 60.0), flexible box\n";
  ang.SetSurveyCoordinates(60.0,0.0);
  tmp_pix.SetResolution(128);
  tmp_pix.SetPixnumFromAng(ang);
  for (uint16_t resolution=128;
       resolution>=Stomp::HPixResolution;resolution/=2) {
    tmp_pix.SetToSuperPix(resolution);

    uint32_t y_min, y_max, n_pixel;
    std::vector<uint32_t> x_min, x_max;
    tmp_pix.XYBounds(theta,x_min,x_max,y_min,y_max,true);
    n_pixel = 0;
    for (uint32_t y=y_min,m=0;y<=y_max;y++,m++)
      n_pixel += x_max[m] - x_min[m] + 1;

    Stomp::PixelVector annulus_pix;
    tmp_pix.WithinRadius(theta,annulus_pix);
    std::cout << "\t" << resolution <<
        ", Checked " << n_pixel << " pixels, kept " <<
        annulus_pix.size() << "\n";
  }
}

void PixelAnnulusIntersectionTests() {
  std::cout << "\n";
  std::cout << "*********************************\n";
  std::cout << "*** Annulus Intersection Test ***\n";
  std::cout << "*********************************\n";

  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::AngularCoordinate pix_ang;
  Stomp::Pixel tmp_pix(ang, 4);
  tmp_pix.Ang(pix_ang);
  ang.SetSurveyCoordinates(pix_ang.Lambda(), tmp_pix.EtaMin());
  double pixel_radius = pix_ang.AngularDistance(ang);
  Stomp::StompWatch stomp_watch;
  std::cout << "Pixel at " << tmp_pix.Resolution() <<
    " resolution; pixel radius = " << pixel_radius << "\n";

  // First, we try an annulus that should fully enclose the pixel.
  tmp_pix.Ang(ang);
  double theta = 10.0*pixel_radius;
  std::cout << "Annulus containing pixel: " << 0.0 << " - " << theta <<
    "\n";
  if (tmp_pix.IntersectsAnnulus(ang, 0.0, theta) == 1) {
    std::cout << "\tGood, a circle contains the pixel.\n";
  } else {
    std::cout << "\tBad, a circle doesn't contain the pixel.\n";
    std::cout << "\tIntersectsAnnulus = " <<
      static_cast<int>(tmp_pix.IntersectsAnnulus(ang, 0.0, theta)) <<
      "; should == 1\n";
    if (tmp_pix.Contains(ang)) {
      std::cout << "\t\tContained in pixel.\n";
    } else {
      std::cout << "\t\tOutside pixel.\n";
    }
    std::cout << "\t\tEdge distances: " <<
      asin(sqrt(tmp_pix.NearEdgeDistance(ang)))*Stomp::RadToDeg <<
      " - " <<
      asin(sqrt(tmp_pix.FarEdgeDistance(ang)))*Stomp::RadToDeg << "\n";
    std::cout << "\t\tCorner distances: " <<
      asin(sqrt(tmp_pix.NearEdgeDistance(ang)))*Stomp::RadToDeg <<
      " - " <<
      asin(sqrt(tmp_pix.FarEdgeDistance(ang)))*Stomp::RadToDeg << "\n";
  }
  stomp_watch.StartTimer();
  int8_t intersects_annulus = tmp_pix.IntersectsAnnulus(ang, 0.0, theta);
  stomp_watch.StopTimer();
  std::cout <<
    "\t\tElapsed Time: " << stomp_watch.ElapsedTime() << " seconds.\n";

  // Now, we shrink the radius, but keep the center at the pixel center.  This
  // should indicate that the pixel intersects, but doesn't contain (== -1).
  tmp_pix.Ang(ang);
  theta = 0.1*pixel_radius;
  std::cout << "Circle inside pixel: " << theta << "\n";
  if (tmp_pix.IntersectsAnnulus(ang, 0.0, theta) == -1) {
    std::cout << "\tGood, circle comes back as intersecting.\n";
  } else {
    std::cout << "\tBad, circle doesn't come back as intersecting.\n";
    std::cout << "\tIntersectsAnnulus = " <<
      tmp_pix.IntersectsAnnulus(ang, 0.0, theta) << "; should == -1\n";
    if (tmp_pix.Contains(ang)) {
      std::cout << "\t\tContained in pixel.\n";
    } else {
      std::cout << "\t\tOutside pixel.\n";
    }
    std::cout << "\t\tEdge distances: " <<
      asin(sqrt(tmp_pix.NearEdgeDistance(ang)))*Stomp::RadToDeg <<
      " - " <<
      asin(sqrt(tmp_pix.FarEdgeDistance(ang)))*Stomp::RadToDeg << "\n";
    std::cout << "\t\tCorner distances: " <<
      asin(sqrt(tmp_pix.NearEdgeDistance(ang)))*Stomp::RadToDeg <<
      " - " <<
      asin(sqrt(tmp_pix.FarEdgeDistance(ang)))*Stomp::RadToDeg << "\n";
  }
  stomp_watch.StartTimer();
  intersects_annulus = tmp_pix.IntersectsAnnulus(ang, 0.0, theta);
  stomp_watch.StopTimer();
  std::cout <<
    "\t\tElapsed Time: " << stomp_watch.ElapsedTime() << " seconds.\n";

  // Inscribed test again, but we make it an annulus.
  double theta_min = 0.8*pixel_radius;
  double theta_max = 1.25*pixel_radius;
  std::cout << "Annulus inside pixel: " << theta_min << " - " <<
    theta_max << "\n";
  if (tmp_pix.IntersectsAnnulus(ang, theta_min, theta_max) == -1) {
    std::cout << "\tGood, annulus comes back as intersecting.\n";
  } else {
    std::cout << "\tBad, annulus doesn't come back as intersecting.\n";
    std::cout << "\tIntersectsAnnulus = " <<
      tmp_pix.IntersectsAnnulus(ang, theta_min, theta_max) <<
      "; should == -1\n";
    if (tmp_pix.Contains(ang)) {
      std::cout << "\t\tContained in pixel.\n";
    } else {
      std::cout << "\t\tOutside pixel.\n";
    }
    std::cout << "\t\tEdge distances: " <<
      asin(sqrt(tmp_pix.NearEdgeDistance(ang)))*Stomp::RadToDeg <<
      " - " <<
      asin(sqrt(tmp_pix.FarEdgeDistance(ang)))*Stomp::RadToDeg << "\n";
    std::cout << "\t\tCorner distances: " <<
      asin(sqrt(tmp_pix.NearEdgeDistance(ang)))*Stomp::RadToDeg <<
      " - " <<
      asin(sqrt(tmp_pix.FarEdgeDistance(ang)))*Stomp::RadToDeg << "\n";
  }
  stomp_watch.StartTimer();
  intersects_annulus = tmp_pix.IntersectsAnnulus(ang, theta_min, theta_max);
  stomp_watch.StopTimer();
  std::cout <<
    "\t\tElapsed Time: " << stomp_watch.ElapsedTime() << " seconds.\n";

  // Circle outside the pixel.
  ang.SetSurveyCoordinates(tmp_pix.Lambda(), tmp_pix.Eta()+10.0*pixel_radius);
  theta = 0.1*pixel_radius;
  std::cout << "Circle " << pix_ang.AngularDistance(ang) <<
    " from pixel, " << theta << " radius\n";
  if (tmp_pix.IntersectsAnnulus(ang, 0.0, theta) == 0) {
    std::cout << "\tGood, circle comes back as outside.\n";
  } else {
    std::cout << "\tBad, circle doesn't come back as outside.\n";
    std::cout << "\tIntersectsAnnulus = " <<
      tmp_pix.IntersectsAnnulus(ang, 0.0, theta) << "; should == 0\n";
    if (tmp_pix.Contains(ang)) {
      std::cout << "\t\tContained in pixel.\n";
    } else {
      std::cout << "\t\tOutside pixel.\n";
    }
    std::cout << "\t\tEdge distances: " <<
      asin(sqrt(tmp_pix.NearEdgeDistance(ang)))*Stomp::RadToDeg <<
      " - " <<
      asin(sqrt(tmp_pix.FarEdgeDistance(ang)))*Stomp::RadToDeg << "\n";
    std::cout << "\t\tCorner distances: " <<
      asin(sqrt(tmp_pix.NearEdgeDistance(ang)))*Stomp::RadToDeg <<
      " - " <<
      asin(sqrt(tmp_pix.FarEdgeDistance(ang)))*Stomp::RadToDeg << "\n";
  }
  stomp_watch.StartTimer();
  intersects_annulus = tmp_pix.IntersectsAnnulus(ang, 0.0, theta);
  stomp_watch.StopTimer();
  std::cout <<
    "\t\tElapsed Time: " << stomp_watch.ElapsedTime() << " seconds.\n";

  // Annulus outside the pixel.
  ang.SetSurveyCoordinates(tmp_pix.Lambda(), tmp_pix.Eta()+10.0*pixel_radius);
  theta_min = 0.1*pixel_radius;
  theta_max = 0.2*pixel_radius;
  std::cout << "Annulus " << pix_ang.AngularDistance(ang) <<
    " from pixel, " << theta_min << " - " << theta_max << "\n";
  if (tmp_pix.IntersectsAnnulus(ang, theta_min, theta_max) == 0) {
    std::cout << "\tGood, annulus comes back as outside.\n";
  } else {
    std::cout << "\tBad, annulus doesn't come back as outside.\n";
    std::cout << "\tIntersectsAnnulus = " <<
      tmp_pix.IntersectsAnnulus(ang, theta_min, theta_max) << "; should == 0\n";
    if (tmp_pix.Contains(ang)) {
      std::cout << "\t\tContained in pixel.\n";
    } else {
      std::cout << "\t\tOutside pixel.\n";
    }
    std::cout << "\t\tEdge distances: " <<
      asin(sqrt(tmp_pix.NearEdgeDistance(ang)))*Stomp::RadToDeg <<
      " - " <<
      asin(sqrt(tmp_pix.FarEdgeDistance(ang)))*Stomp::RadToDeg << "\n";
    std::cout << "\t\tCorner distances: " <<
      asin(sqrt(tmp_pix.NearEdgeDistance(ang)))*Stomp::RadToDeg <<
      " - " <<
      asin(sqrt(tmp_pix.FarEdgeDistance(ang)))*Stomp::RadToDeg << "\n";
  }
  stomp_watch.StartTimer();
  intersects_annulus = tmp_pix.IntersectsAnnulus(ang, theta_min, theta_max);
  stomp_watch.StopTimer();
  std::cout <<
    "\t\tElapsed Time: " << stomp_watch.ElapsedTime() << " seconds.\n";

  // Annulus with center outside the pixel but which should contain the pixel.
  ang.SetSurveyCoordinates(tmp_pix.Lambda(), tmp_pix.Eta()+10.0*pixel_radius);
  theta_min = 0.1*pixel_radius;
  theta_max = 20.0*pixel_radius;
  std::cout << "Annulus " << pix_ang.AngularDistance(ang) <<
    " from pixel, " <<
    theta_min << " - " << theta_max << "\n";
  if (tmp_pix.IntersectsAnnulus(ang, theta_min, theta_max) == 1) {
    std::cout << "\tGood, containing annulus listed as containing.\n";
  } else {
    std::cout << "\tBad, containing annulus listed as not containing.\n";
    std::cout << "\tIntersectsAnnulus = " <<
      tmp_pix.IntersectsAnnulus(ang, theta_min, theta_max) << "; should == 1\n";
    if (tmp_pix.Contains(ang)) {
      std::cout << "\t\tContained in pixel.\n";
    } else {
      std::cout << "\t\tOutside pixel.\n";
    }
    std::cout << "\t\tEdge distances: " <<
      asin(sqrt(tmp_pix.NearEdgeDistance(ang)))*Stomp::RadToDeg <<
      " - " <<
      asin(sqrt(tmp_pix.FarEdgeDistance(ang)))*Stomp::RadToDeg << "\n";
    std::cout << "\t\tCorner distances: " <<
      asin(sqrt(tmp_pix.NearEdgeDistance(ang)))*Stomp::RadToDeg <<
      " - " <<
      asin(sqrt(tmp_pix.FarEdgeDistance(ang)))*Stomp::RadToDeg << "\n";
  }
  stomp_watch.StartTimer();
  intersects_annulus = tmp_pix.IntersectsAnnulus(ang, theta_min, theta_max);
  stomp_watch.StopTimer();
  std::cout <<
    "\t\tElapsed Time: " << stomp_watch.ElapsedTime() << " seconds.\n";

  // Annulus with center outside the pixel but which should intersect the pixel.
  ang.SetSurveyCoordinates(tmp_pix.Lambda(), tmp_pix.Eta()+2.0*pixel_radius);
  theta_min = 1.25*pixel_radius;
  theta_max = 2.5*pixel_radius;
  std::cout << "Annulus " << pix_ang.AngularDistance(ang) <<
    " from pixel, " << theta_min << " - " << theta_max << "\n";
  if (tmp_pix.IntersectsAnnulus(ang, theta_min, theta_max) == -1) {
    std::cout << "\tGood, intersecting annulus listed as intersecting.\n";
  } else {
    std::cout << "\tBad, intersecting annulus is not intersecting.\n";
    std::cout << "\tIntersectsAnnulus = " <<
      tmp_pix.IntersectsAnnulus(ang, theta_min, theta_max) <<
      "; should == -1\n";
    if (tmp_pix.Contains(ang)) {
      std::cout << "\t\tContained in pixel.\n";
    } else {
      std::cout << "\t\tOutside pixel.\n";
    }
    std::cout << "\t\tEdge distances: " <<
      asin(sqrt(tmp_pix.NearEdgeDistance(ang)))*Stomp::RadToDeg <<
      " - " <<
      asin(sqrt(tmp_pix.FarEdgeDistance(ang)))*Stomp::RadToDeg << "\n";
    std::cout << "\t\tCorner distances: " <<
      asin(sqrt(tmp_pix.NearEdgeDistance(ang)))*Stomp::RadToDeg <<
      " - " <<
      asin(sqrt(tmp_pix.FarEdgeDistance(ang)))*Stomp::RadToDeg << "\n";
  }
  stomp_watch.StartTimer();
  intersects_annulus = tmp_pix.IntersectsAnnulus(ang, theta_min, theta_max);
  stomp_watch.StopTimer();
  std::cout <<
    "\t\tElapsed Time: " << stomp_watch.ElapsedTime() << " seconds.\n";
}

void PixelUnitTests(bool run_all_tests) {
  void PixelBasicTests();
  void PixelStripeTests();
  void PixelXYTests();
  void PixelBoundTests();
  void PixelWithinRadiusTests();
  void PixelAnnulusIntersectionTests();

  if (run_all_tests) FLAGS_all_pixel_tests = true;

  // Check that the hierarchical resolution scaling routines work properly.
  if (FLAGS_all_pixel_tests || FLAGS_pixel_basic_tests)
    PixelBasicTests();

  // Check that the routines for generating the stripe from a pixel work.
  if (FLAGS_all_pixel_tests || FLAGS_pixel_stripe_tests) PixelStripeTests();

  // Check that the routines for generating the X-Y bounds of a given
  // region work correctly.
  if (FLAGS_all_pixel_tests || FLAGS_pixel_xy_tests) PixelXYTests();

  // Quick checks for the various routines to return the pixel bounds in the
  // various coordinate systems.
  if (FLAGS_all_pixel_tests || FLAGS_pixel_bound_tests) PixelBoundTests();

  // Check the routines for taking the X-Y bounds and returning a list of pixels
  // with the specified radius.
  if (FLAGS_all_pixel_tests || FLAGS_pixel_within_radius_tests)
    PixelWithinRadiusTests();

  // Check the routines for determining whether or not annuli intersect pixels.
  if (FLAGS_all_pixel_tests || FLAGS_pixel_annulus_intersection_tests)
    PixelAnnulusIntersectionTests();
}
