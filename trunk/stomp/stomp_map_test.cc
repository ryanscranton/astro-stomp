#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <gflags/gflags.h>
#include "stomp_core.h"
#include "stomp_angular_coordinate.h"
#include "stomp_pixel.h"
#include "stomp_geometry.h"
#include "stomp_map.h"

void MapBasicTests() {
  // Ok, now we're ready to start playing with the Stomp::Map interfaces.  We'll
  // start by making a vector of pixels at the same resolution that cover a
  // patch of sky.  Then, we invoke the Stomp::Map constructor, which should
  // resolve that vector into a smaller group of pixels at different
  // resolutions that cover the same area.
  std::cout << "\n";
  std::cout << "***********************\n";
  std::cout << "*** Map Basic Tests ***\n";
  std::cout << "***********************\n";
  double theta = 3.0;
  double true_circle_area =
    (1.0 - cos(theta*Stomp::DegToRad))*
    2.0*Stomp::Pi*Stomp::StradToDeg;
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);

  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);

  uint32_t original_annulus_size = annulus_pix.size();
  std::cout << "\nStomp::Map routines\n";
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);
  std::cout << "\tMade Stomp::Map with " << stomp_map->Size() <<
    " pixels from initial vector of " <<
    original_annulus_size << " pixels.\n";

  // We won't match the true area since we're using relatively coarse pixels,
  // but it should be within a tenth of a square degree or so.
  std::cout << "\tTotal Area: " << stomp_map->Area() << " square degrees.\n";
  std::cout << "\tOriginal Area: " << tmp_pix.Area()*original_annulus_size <<
      " square degrees.\n";
  std::cout << "\tTrue Annulus Area: " << true_circle_area <<
    " square degrees.\n";
}

void MapWriteTests() {
  std::cout << "\n";
  std::cout << "***********************\n";
  std::cout << "*** Map Write Tests ***\n";
  std::cout << "***********************\n";

  double theta = 3.0;
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  std::string output_file_name = "Map.pix";
  std::cout << "\t Writing Stomp::Map to " << output_file_name.c_str() << "\n";
  stomp_map->Write(output_file_name);
}

void MapReadTests() {
  std::cout << "\n";
  std::cout << "**********************\n";
  std::cout << "*** Map Read Tests ***\n";
  std::cout << "**********************\n";

  double theta = 3.0;
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  for (Stomp::PixelIterator iter=annulus_pix.begin();
       iter!=annulus_pix.end();++iter) iter->SetWeight(1.0*iter->Superpixnum());

  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  std::string output_file_name = "StompReadMap.pix";
  std::cout << "\tWriting Stomp::Map to " << output_file_name << "\n";
  stomp_map->Write(output_file_name);
  std::cout << "\t\tDone.\n";

  std::cout << "\tReading Stomp::Map from " << output_file_name << "\n";
  Stomp::Map* read_stomp_map = new Stomp::Map(output_file_name);
  std::cout << "\t\tDone.\n";

  std::cout << "\nChecking Map parameters...\n";
  std::cout << "\tArea: " << read_stomp_map->Area() << " (" <<
    stomp_map->Area() << ")\n";
  std::cout << "\tSize: " << read_stomp_map->Size() << " (" <<
    stomp_map->Size() << ")\n";
  std::cout << "\tWeight: " << read_stomp_map->MinWeight() << " - " <<
    read_stomp_map->MaxWeight() << " (" << stomp_map->MinWeight() << " - " <<
    stomp_map->MaxWeight() << ")\n";
  std::cout << "\tResolution: " << read_stomp_map->MinResolution() << " - " <<
    read_stomp_map->MaxResolution() << " (" << stomp_map->MinResolution() <<
    " - " << stomp_map->MaxResolution() << ")\n";
  std::cout << "\tResolution break-down:\n";
  for (uint32_t resolution=Stomp::HPixResolution, i=0;
       i<Stomp::ResolutionLevels;resolution*=2, i++)
    std::cout << "\t\t" << resolution << ": " <<
      read_stomp_map->PixelCount(resolution) << " (" <<
      stomp_map->PixelCount(resolution) << ")\n";
}

void MapPixelizationTests() {
  // Now we try generating a Map from a GeometricBound object
  std::cout << "\n";
  std::cout << "******************************\n";
  std::cout << "*** Map Pixelization Tests ***\n";
  std::cout << "******************************\n";

  Stomp::AngularCoordinate ang(20.0, 0.0, Stomp::AngularCoordinate::Survey);
  double radius = 3.0;
  Stomp::CircleBound circle(ang, radius);

  Stomp::Map* stomp_map = new Stomp::Map();

  std::cout << "Pixelizing at full resolution...\n";
  uint8_t max_resolution_level = Stomp::MaxPixelLevel;
  uint8_t resolution_level =
    stomp_map->_FindStartingResolutionLevel(circle.Area());
  std::cout << "\tPixelization resolution level range: " <<
    static_cast<int>(resolution_level) << " - " <<
    static_cast<int>(max_resolution_level) << "\n";

  if (stomp_map->PixelizeBound(circle, 1.0)) {
    std::cout << "\tPixelization success!\n\t\tArea Check: Real: " <<
      circle.Area() << ", Pixelized: " << stomp_map->Area() << "; " <<
      stomp_map->Size() << " pixels.\n";
  } else {
    std::cout << "Pixelization failed...\n";
  }

  uint32_t resolution = 2048;
  std::cout << "\nPixelizing at coarser resolution level: " <<
    static_cast<int>(resolution_level) << " - " <<
    static_cast<int>(Stomp::MostSignificantBit(resolution)) << "\n";
  if (stomp_map->PixelizeBound(circle, 1.0, resolution)) {
    std::cout << "\tPixelization success!\n\t\tArea Check: Real: " <<
      circle.Area() << ", Pixelized: " << stomp_map->Area() << "; " <<
      stomp_map->Size() << " pixels.\n";
  } else {
    std::cout << "Pixelization failed...\n";
  }
}

void MapCoverTests() {
  std::cout << "\n";
  std::cout << "***********************\n";
  std::cout << "*** Map Cover Tests ***\n";
  std::cout << "***********************\n";
  double theta = 3.0;
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  Stomp::PixelVector superpix;
  stomp_map->Coverage(superpix);
  std::cout << "\tStomp::Map covers " << superpix.size() << " superpixels:\n";
  for (Stomp::PixelIterator iter=superpix.begin();iter!=superpix.end();++iter) {
    uint32_t k = iter->Superpixnum();
    std::cout << "\t\t" << k << ": " << stomp_map->Area(k) <<
      " sq. degrees, " << stomp_map->MinResolution(k) <<
      " <= resolution <= " << stomp_map->MaxResolution(k) << "\n";
  }

  Stomp::Map tmp_map;
  bool made_covering = stomp_map->Covering(tmp_map, 10*superpix.size());
  if (made_covering) {
    std::cout << "\tMade covering Stomp::Map with " << tmp_map.Size() <<
      " (" << 10*superpix.size() << ") pixels\n";
  } else {
    std::cout << "Failed to make covering Stomp::Map with <" <<
      10*superpix.size() << " pixels (" << tmp_map.Size() << ")\n";
  }
}

void MapIteratorTests() {
  std::cout << "\n";
  std::cout << "**************************\n";
  std::cout << "*** Map Iterator Tests ***\n";
  std::cout << "**************************\n";
  double theta = 3.0;
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  for (Stomp::PixelIterator iter=annulus_pix.begin();
       iter!=annulus_pix.end();++iter) iter->SetWeight(1.0*iter->Superpixnum());
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  std::cout << "Initial Map parameters:\n";
  Stomp::PixelVector superpix;
  stomp_map->Coverage(superpix);
  std::cout << "\t" << superpix.size() << " superpixels.\n";
  for (Stomp::PixelIterator iter=superpix.begin();
       iter!=superpix.end();++iter)
    std::cout << "\t\t" << iter->Superpixnum() << ": " <<
      stomp_map->Size(iter->Superpixnum()) << " pixels, " <<
      stomp_map->Area(iter->Superpixnum()) << " sq. degrees.\n";

  uint32_t counter = 0;
  uint32_t size = 0;
  double area = 0.0;
  double weight_min = 1.0e30;
  double weight_max = -1.0e30;
  uint32_t resolution_min = Stomp::MaxPixelResolution;
  uint32_t resolution_max = Stomp::HPixResolution;
  Stomp::ResolutionDict pixel_count;
  std::cout << "\nAttempting to iterate through " << stomp_map->Size() <<
    " pixels...\n";
  for (Stomp::MapIterator iter=stomp_map->Begin();
       iter!=stomp_map->End();stomp_map->Iterate(&iter), counter++) {
    uint32_t resolution = iter.second->Resolution();
    uint32_t resolution_counter = 0;
    for (uint32_t resolution_mask=Stomp::HPixResolution, i=0;
	 i<Stomp::ResolutionLevels;resolution_mask*=2, i++)
      if (resolution & resolution_mask) resolution_counter++;
    if (resolution_counter == 1) {
      size++;
      area += iter.second->Area();
      double weight = iter.second->Weight();
      if (weight < weight_min) weight_min = weight;
      if (weight > weight_max) weight_max = weight;
      if (resolution < resolution_min) resolution_min = resolution;
      if (resolution > resolution_max) resolution_max = resolution;
      pixel_count[iter.second->Resolution()] += 1;
    } else {
      std::cout << "\t" << counter << " Bogus pixel: resolution = " <<
	resolution << ", weight = " << iter.second->Weight() <<
	", Superpixnum = " << iter.first << "\n";
    }
  }

  std::cout << "\nChecking Map parameters...\n";
  std::cout << "\tArea: " << area << " (" << stomp_map->Area() << ")\n";
  std::cout << "\tSize: " << size << " (" << stomp_map->Size() << ")\n";
  std::cout << "\tWeight: " << weight_min << " - " << weight_max <<
    " (" << stomp_map->MinWeight() << " - " << stomp_map->MaxWeight() << ")\n";
  std::cout << "\tResolution: " << resolution_min << " - " << resolution_max <<
    " (" << stomp_map->MinResolution() << " - " <<
    stomp_map->MaxResolution() << ")\n";
  std::cout << "\tResolution break-down:\n";
  for (uint32_t resolution=Stomp::HPixResolution, i=0;
       i<Stomp::ResolutionLevels;resolution*=2, i++)
    std::cout << "\t\t" << resolution << ": " << pixel_count[resolution] <<
      " (" << stomp_map->PixelCount(resolution) << ")\n";
}

void MapLocationTests() {
  std::cout << "\n";
  std::cout << "**************************\n";
  std::cout << "*** Map Location Tests ***\n";
  std::cout << "**************************\n";
  double theta = 3.0;
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  // Ok, now we'll check some points that may or may not be in the map.
  ang.SetSurveyCoordinates(60.0, 0.0);
  double weight = 0.0;

  if (stomp_map->FindLocation(ang, weight)) {
    std::cout <<
      "\tGood. Found point at (60,0) which was the center of the map\n";
    std::cout << "\tThe weight here is " << weight << " (1.0).\n";
  } else {
    std::cout << "\tNot good. That point (60,0) was the center of the map\n";
  }

  ang.SetSurveyCoordinates(0.0, 0.0);
  if (stomp_map->FindLocation(ang, weight)) {
    std::cout << "\tNot good. That point (0,0) was well outside the map\n";
  } else {
    std::cout << "\tGood. That point (0,0) was well outside the map\n";
    std::cout << "\t\tThe weight here is " << weight << " (1.0).\n";
  }

  ang.SetSurveyCoordinates(59.0, 1.0);
  if (stomp_map->FindLocation(ang, weight)) {
    std::cout << "\tGood. Found point at (59,1) which was within the map\n";
    std::cout << "\t\tThe weight here is " << weight << " (1.0).\n";
  } else {
    std::cout << "\tNot good. That point (59,1) was within the map\n";
  }

  ang.SetSurveyCoordinates(60.0, 10.0);
  if (stomp_map->FindLocation(ang, weight)) {
    std::cout << "\tNot good. That point (60,10) was outside the map\n";
  } else {
    std::cout << "\tGood. That point (60,10) was well outside the map\n";
    std::cout << "\t\tThe weight here is " << weight << " (1.0).\n";
  }
}

void MapUnmaskedFractionTests() {
  // Ok, now to test the unmasked fraction code.  We'll scroll through the
  // the superpixels and make sure that the area returned matches the area
  // we've already established for them.
  std::cout << "\n";
  std::cout << "***********************************\n";
  std::cout << "*** Map Unmasked Fraction Tests ***\n";
  std::cout << "***********************************\n";
  double theta = 3.0;
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  Stomp::PixelVector superpix;
  stomp_map->Coverage(superpix);
  for (Stomp::PixelIterator iter=superpix.begin();iter!=superpix.end();++iter) {
    double unmasked_fraction = stomp_map->FindUnmaskedFraction(*iter);
    std::cout << "\t" << iter->Superpixnum() << ": " <<
      unmasked_fraction*Stomp::HPixArea << " (" <<
      stomp_map->Area(iter->Superpixnum()) << ") square degrees.\n";
  }
}

void MapContainsTests() {
  // Ok, now we want to test the various routines for determining whether or
  // not a map is contained within another.
  std::cout << "\n";
  std::cout << "*************************\n";
  std::cout << "*** Map Contain Tests ***\n";
  std::cout << "*************************\n";
  double theta = 3.0;
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  tmp_pix.WithinRadius(theta/5.0, annulus_pix);
  Stomp::Map* stomp_map_inside = new Stomp::Map(annulus_pix);

  ang.SetSurveyCoordinates(0.0,0.0);
  tmp_pix.SetResolution(256);
  tmp_pix.SetPixnumFromAng(ang);
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map_outside = new Stomp::Map(annulus_pix);

  ang.SetSurveyCoordinates(60.0,0.5*theta);
  tmp_pix.SetResolution(256);
  tmp_pix.SetPixnumFromAng(ang);
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map_partial = new Stomp::Map(annulus_pix);

  if (stomp_map->Contains(*stomp_map_inside)) {
    std::cout << "Good.  Map inside our Map comes back as contained (" <<
      static_cast<int>(stomp_map->FindUnmaskedStatus(*stomp_map_inside)) <<
      ").\n";
  } else {
    std::cout << "Bad.  Map inside our Map comes back as not contained (" <<
      static_cast<int>(stomp_map->FindUnmaskedStatus(*stomp_map_inside)) <<
      ").\n";
  }

  if (!stomp_map->Contains(*stomp_map_outside)) {
    std::cout << "Good.  Map outside our Map comes back as not contained (" <<
      static_cast<int>(stomp_map->FindUnmaskedStatus(*stomp_map_outside)) <<
      ").\n";
  } else {
    std::cout << "Bad.  Map outside our Map comes back as contained (" <<
      static_cast<int>(stomp_map->FindUnmaskedStatus(*stomp_map_outside)) <<
      ").\n";
  }

  if (!stomp_map->Contains(*stomp_map_partial)) {
    std::cout <<
      "Good.  Map partially outside our Map comes back as not contained (" <<
      static_cast<int>(stomp_map->FindUnmaskedStatus(*stomp_map_partial)) <<
      ").\n";
  } else {
    std::cout <<
      "Bad.  Map partially outside our Map comes back as contained (" <<
      static_cast<int>(stomp_map->FindUnmaskedStatus(*stomp_map_partial)) <<
      ").\n";
  }
}

void MapRandomPointsTests() {
  // Alright, now we check the random position generator.  This should give
  // us back a fixed number of randomly selected positions within our original
  // map, so we'll check that that's happening.
  std::cout << "\n";
  std::cout << "*******************************\n";
  std::cout << "*** Map Random Points Tests ***\n";
  std::cout << "*******************************\n";
  double theta = 3.0;
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  uint32_t n_random = 10000;
  std::cout << "\tGenerating " << n_random << " points...\n";
  Stomp::AngularVector rand_ang;
  stomp_map->GenerateRandomPoints(rand_ang,n_random);

  if (rand_ang.size() != n_random) {
    std::cout << "\t\tRandom point array size doesn't match requested size!\n";
    exit(1);
  }

  uint32_t n_found = 0;
  double weight = 0.0;
  for (Stomp::AngularIterator iter=rand_ang.begin();iter!=rand_ang.end();++iter)
    if (stomp_map->FindLocation(*iter,weight)) n_found++;
  std::cout << "\tVerified that " << n_found << "/" << n_random <<
      " points within map.\n";

  ang.SetSurveyCoordinates(62.0,2.0);
  tmp_pix.SetResolution(256);
  tmp_pix.SetPixnumFromAng(ang);
  tmp_pix.WithinRadius(theta,annulus_pix);
  Stomp::Map* stomp_map_nearby = new Stomp::Map(annulus_pix);

  ang.SetSurveyCoordinates(0.0,0.0);
  tmp_pix.SetResolution(256);
  tmp_pix.SetPixnumFromAng(ang);
  tmp_pix.WithinRadius(theta,annulus_pix);
  Stomp::Map* stomp_map_faraway = new Stomp::Map(annulus_pix);

  // Now we'll check the random points from our original map against the two
  // new maps and see how much overlap we've got.
  std::cout << "\tTesting against random positions from original map:\n";
  n_found = 0;
  for (Stomp::AngularIterator iter=rand_ang.begin();iter!=rand_ang.end();++iter)
    if (stomp_map_nearby->FindLocation(*iter,weight)) n_found++;
  std::cout << "\t\tFound " << n_found << "/" << n_random <<
      " points within nearby map.\n";

  n_found = 0;
  for (Stomp::AngularIterator iter=rand_ang.begin();iter!=rand_ang.end();++iter)
    if (stomp_map_faraway->FindLocation(*iter,weight)) n_found++;
  std::cout << "\t\tFound " << n_found << "/" << n_random <<
      " points within far away map.\n";
}

void MapMultiMapTests() {
  // Ok, now we want to test the various routines for working with multiple
  // stomp maps.  We'll use the same basic routines to generate two new maps,
  // one that should partially overlap our original map and one that is well
  // away from the other two.
  std::cout << "\n";
  std::cout << "***************************\n";
  std::cout << "*** Map Multi-Map Tests ***\n";
  std::cout << "***************************\n";
  uint32_t pixelate_resolution = 32768;
  double radius = 3.0;
  std::cout << "\tCenter circle:\n\t";
  Stomp::AngularCoordinate center_ang(20.0, 0.0,
				      Stomp::AngularCoordinate::Survey);
  Stomp::CircleBound center_circle(center_ang, radius);
  Stomp::Map* stomp_map =
    new Stomp::Map(center_circle, 1.0, pixelate_resolution, true);
  std::string output_file = "CenterCircle.map";
  stomp_map->Write(output_file);

  std::cout << "\tNearby circle:\n\t";
  Stomp::AngularCoordinate nearby_ang(22.0, 2.0,
				      Stomp::AngularCoordinate::Survey);
  Stomp::CircleBound nearby_circle(nearby_ang, radius);
  Stomp::Map* nearby_map =
    new Stomp::Map(nearby_circle, 1.0, pixelate_resolution, true);
  output_file = "NearbyCircle.map";
  nearby_map->Write(output_file);

  std::cout << "\tFaraway circle:\n\t";
  Stomp::AngularCoordinate faraway_ang(0.0, 0.0,
				      Stomp::AngularCoordinate::Survey);
  Stomp::CircleBound faraway_circle(faraway_ang, radius);
  Stomp::Map* faraway_map =
    new Stomp::Map(faraway_circle, 1.0, pixelate_resolution, true);
  output_file = "FarawayCircle.map";
  faraway_map->Write(output_file);

  // Now, we'll test the routines for doing logical AND, OR and NOT with the
  // various maps.
  std::cout << "\n\t**********************\n";
  std::cout << "\t*** Ingestion Test ***\n";
  std::cout << "\t**********************\n";
  Stomp::Map* tmp_map = new Stomp::Map(center_circle, 1.0, pixelate_resolution);
  tmp_map->IngestMap(*nearby_map, false);
  std::cout << "\t\tOriginal map area: " << stomp_map->Area() <<
    " sq. degrees.\n";
  std::cout << "\t\tNearby: Adding " << nearby_map->Area() <<
      " sq. degree map to original map\n";
  std::cout << "\t\t\tNew Map: " << tmp_map->Area() << " sq. degrees.\n";
  output_file = "NearbyIngestion.map";
  tmp_map->Write(output_file);
  delete tmp_map;

  tmp_map = new Stomp::Map(center_circle, 1.0, pixelate_resolution);
  tmp_map->IngestMap(*faraway_map, false);
  std::cout << "\t\tFar Away: Adding " << faraway_map->Area() <<
      " sq. degree map to original map.\n";
  std::cout << "\t\t\tNew Map: " << tmp_map->Area() << " (" <<
      stomp_map->Area()+faraway_map->Area() << ") sq. degrees.\n";
  output_file = "FarawayIngestion.map";
  tmp_map->Write(output_file);
  delete tmp_map;

  std::cout << "\n\t*************************\n";
  std::cout << "\t*** Intersection Test ***\n";
  std::cout << "\t*************************\n";
  tmp_map = new Stomp::Map(center_circle, 1.0, pixelate_resolution);
  if (tmp_map->IntersectMap(*nearby_map)) {
    std::cout << "\t\tNearby intersection area: " <<
        tmp_map->Area() << " sq. degrees.\n";
    output_file = "NearbyIntersection.map";
    tmp_map->Write(output_file);
  } else {
    std::cout << "\t\tThis is bad," <<
        " there should have been some intersecting area.\n";
  }
  delete tmp_map;

  tmp_map = new Stomp::Map(center_circle, 1.0, pixelate_resolution);
  if (tmp_map->IntersectMap(*faraway_map)) {
    std::cout << "\t\tBad. Far away intersection area: " <<
        tmp_map->Area() << " sq. degrees.  Should be 0 and leave unchanged.\n";
  } else {
    std::cout <<
      "\t\tGood, no intersecting area as expected for Far Away map.\n";
  }
  delete tmp_map;

  std::cout << "\n\t**********************\n";
  std::cout << "\t*** Exclusion Test ***\n";
  std::cout << "\t**********************\n";
  tmp_map = new Stomp::Map(center_circle, 1.0, pixelate_resolution);
  if (tmp_map->ExcludeMap(*nearby_map, false)) {
    std::cout << "\t\tArea remaining after excluding nearby map: " <<
        tmp_map->Area() << " sq. degrees.\n";
  } else {
    std::cout << "\t\tThis is bad," <<
        " there should have been some area left over.\n";
  }
  output_file = "NearbyExclusion.map";
  tmp_map->Write(output_file);
  delete tmp_map;

  tmp_map = new Stomp::Map(center_circle, 1.0, pixelate_resolution);
  if (tmp_map->ExcludeMap(*faraway_map, false)) {
    std::cout << "\t\tArea remaining after excluding faraway map: " <<
        tmp_map->Area() << " sq. degrees.\n";
  } else {
    std::cout <<
      "\t\tBad, we should have recovered all of the original area.\n";
  }
  output_file = "FarawayExclusion.map";
  tmp_map->Write(output_file);
  delete tmp_map;

  std::cout << "\n\t*********************\n";
  std::cout << "\t*** Addition Test ***\n";
  std::cout << "\t*********************\n";
  tmp_map = new Stomp::Map(center_circle, 3.0, pixelate_resolution);
  if (tmp_map->AddMap(*nearby_map, true)) {
    std::cout << "\t\tNearby overlapping area: " <<
        tmp_map->Area() << " sq. degrees.\n";
    std::cout << "\t\tAverage weight: " << tmp_map->AverageWeight() <<
        " (" << tmp_map->MinWeight() << " - " << tmp_map->MaxWeight() << ")\n";
  } else {
    std::cout << "\t\tThis is bad," <<
        " there should have been some area left over.\n";
  }
  output_file = "NearbyAdditionDropSingle.map";
  tmp_map->Write(output_file);
  delete tmp_map;

  tmp_map = new Stomp::Map(center_circle, 3.0, pixelate_resolution);
  if (tmp_map->AddMap(*nearby_map, false)) {
    std::cout << "\t\tNearby combined area: " <<
        tmp_map->Area() << " sq. degrees.\n";
    std::cout << "\t\tAverage weight: " << tmp_map->AverageWeight() <<
        " (" << tmp_map->MinWeight() << " - " << tmp_map->MaxWeight() << ")\n";
  } else {
    std::cout << "\t\tThis is bad," <<
        " there should have been some area left over.\n";
  }
  output_file = "NearbyAdditionKeepSingle.map";
  tmp_map->Write(output_file);
  delete tmp_map;

  tmp_map = new Stomp::Map(center_circle, 2.0, pixelate_resolution);
  if (tmp_map->AddMap(*faraway_map, false)) {
    std::cout << "\t\tFaraway combined area: " << tmp_map->Area() <<
        " sq. degrees.\n";
    std::cout << "\t\tAverage weight: " << tmp_map->AverageWeight() <<
        " (" << tmp_map->MinWeight() << " - " << tmp_map->MaxWeight() << ")\n";
  } else {
    std::cout <<
      "\t\tBad, we should have recovered all of the original area.\n";
  }
  output_file = "FarawayKeepSingle.map";
  tmp_map->Write(output_file);
  delete tmp_map;

  std::cout << "\n\t***************************\n";
  std::cout << "\t*** Multiplication Test ***\n";
  std::cout << "\t***************************\n";

  tmp_map = new Stomp::Map(center_circle, 3.0, pixelate_resolution);
  if (tmp_map->MultiplyMap(*nearby_map, true)) {
    std::cout << "\t\tNearby overlapping area: " <<
        tmp_map->Area() << " sq. degrees.\n";
    std::cout << "\t\tAverage weight: " << tmp_map->AverageWeight() << "\n";
  } else {
    std::cout << "\t\tThis is bad," <<
        " there should have been some area left over.\n";
  }
  output_file = "NearbyMultiplicationDropSingle.map";
  tmp_map->Write(output_file);
  delete tmp_map;

  tmp_map = new Stomp::Map(center_circle, 3.0, pixelate_resolution);
  if (tmp_map->MultiplyMap(*nearby_map, false)) {
    std::cout << "\t\tNearby combined area: " <<
        tmp_map->Area() << " sq. degrees.\n";
    std::cout << "\t\tAverage weight: " << tmp_map->AverageWeight() << "\n";
  } else {
    std::cout << "\t\tThis is bad," <<
        " there should have been some area left over.\n";
  }
  output_file = "NearbyMultiplicationKeepSingle.map";
  tmp_map->Write(output_file);
  delete tmp_map;

  tmp_map = new Stomp::Map(center_circle, 2.0, pixelate_resolution);
  if (tmp_map->MultiplyMap(*faraway_map, false)) {
    std::cout << "\t\tFaraway combined area: " << tmp_map->Area() <<
        " sq. degrees.\n";
    std::cout << "\t\tAverage weight: " << tmp_map->AverageWeight() << "\n";
  } else {
    std::cout <<
      "\t\tBad, we should have recovered all of the original area.\n";
  }
  output_file = "FarawayMultiplicationKeepSingle.map";
  tmp_map->Write(output_file);
  delete tmp_map;

  delete nearby_map;
  delete faraway_map;
}

void MapRegionTests() {
  // Ok, now we want to test the various routines for breaking our map into
  // nearly equal sub-regions.
  std::cout << "\n";
  std::cout << "************************\n";
  std::cout << "*** Map Region Tests ***\n";
  std::cout << "************************\n\n";
  double theta = 3.0;
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, 256);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);
  uint16_t n_regions = 10;

  std::cout << "Trying to regionate the map into " << n_regions <<
    " pieces...\n";
  stomp_map->InitializeRegions(n_regions, Stomp::HPixResolution);

  // A better result should be found as we increase the resolution for the
  // region map.
  std::cout << "Now doing it with finer resolution...\n";
  stomp_map->InitializeRegions(n_regions);

  // Now we check to see if which region contains the point used to generate
  // our annulus.
  int16_t center_region = stomp_map->FindRegion(ang);
  std::cout << "Center point is in region " << center_region << "\n";

  // Finally, we want to check our routines for returning maps that include
  // just a given region or all of the map except for a given region.
  Stomp::Map region_map;
  stomp_map->RegionOnlyMap(center_region, region_map);

  std::cout << "Total map area: " << stomp_map->Area() << "\n";
  std::cout << "Region " << center_region << " area: " <<
    region_map.Area() << "\n";
  if (region_map.Contains(ang)) {
    std::cout << "\tAnd the Stomp::Map from that region contains " <<
      "the center point.  Good.\n";
  } else {
    std::cout << "\tThis is bad; the center point should be in this map.\n";
  }

  Stomp::Map excluded_region_map;
  stomp_map->RegionExcludedMap(center_region, excluded_region_map);
  std::cout << "Region " << center_region << " excluded area: " <<
    excluded_region_map.Area() << "\n";
  if (excluded_region_map.Contains(ang)) {
    std::cout << "\tThis is bad; the center point shouldn't be in this map.\n";
  } else {
    std::cout <<
      "\tAnd the Stomp::Map excluding that region doesn't contain " <<
      "the center point.  Good.\n";
  }
}

void MapRegionBoundTests() {
  // Ok, now we want to test the various routines for breaking our map into
  // nearly equal sub-regions where we have user-supplied RegionBounds
  std::cout << "\n";
  std::cout << "*****************************\n";
  std::cout << "*** Map RegionBound Tests ***\n";
  std::cout << "*****************************\n\n";
  Stomp::LatLonBound map_bound(-5.0, 5.0, -5.0, 5.0,
			       Stomp::AngularCoordinate::Equatorial);
  std::cout << "Map bounds: " <<
    map_bound.LatitudeMin() << " - " << map_bound.LatitudeMax() << ", " <<
    map_bound.LongitudeMin() << " - " << map_bound.LongitudeMax() << "\n";

  Stomp::Map* stomp_map = new Stomp::Map(map_bound, 1.0, 256, true);
  uint16_t n_regions = 10;

  // First, we regionate the area using the standard method.
  std::cout << "Trying to regionate the map into " << n_regions <<
    " pieces...\n";
  uint16_t n_regions_final = stomp_map->InitializeRegions(n_regions);
  std::cout << "\tRegionated into " << n_regions_final << " (" <<
    n_regions << ") regions at " << stomp_map->RegionResolution() << "\n";
  for (uint16_t region_iter=0;region_iter<n_regions_final;region_iter++)
    std::cout << "\t\t" << region_iter << ": " <<
      stomp_map->RegionArea(region_iter) << " sq. degrees.\n";

  // Now, we create a RegionBound in the center of the Map and re-regionate
  // with it included.
  Stomp::LatLonBound* single_bound =
    new Stomp::LatLonBound(-1.0, 1.0, -1.0, 1.0,
			   Stomp::AngularCoordinate::Equatorial);
  Stomp::RegionBound region_bound(single_bound);

  std::cout << "Checking RegionBound behavior..\n";
  Stomp::PixelVector coverage_pix;
  stomp_map->Coverage(coverage_pix);
  uint32_t bound_contained = 0;
  uint32_t region_contained = 0;
  for (Stomp::PixelIterator iter=coverage_pix.begin();
       iter!=coverage_pix.end();++iter) {
    if (single_bound->CheckPixel(*iter)) bound_contained++;
    if (region_bound.CheckPixel(*iter)) region_contained++;
  }
  std::cout << "\t" << coverage_pix.size() <<
    " coverage pixels at resolution " << Stomp::HPixResolution << "\n\t\t" <<
    bound_contained << " contained in GeometricBound\n\t\t" <<
    region_contained << " contained in RegionBound\n";
  std::cout << "\tGeometricBound area: " << single_bound->Area() <<
    " sq. deg.\n";
  std::cout << "\tRegionBound area: " << region_bound.BoundArea() <<
    " sq. deg.\n";

  Stomp::RegionBoundVector regionVec;
  regionVec.push_back(region_bound);

  std::cout << "\nRe-Regionating the map with central RegionBound...\n";
  n_regions_final = stomp_map->InitializeRegions(regionVec, n_regions);
  std::cout << "\tRegionated into " << n_regions_final << " (" <<
    n_regions << ") regions at " << stomp_map->RegionResolution() << "\n";
  for (uint16_t region_iter=0;region_iter<n_regions_final;region_iter++)
    std::cout << "\t\t" << region_iter << ": " <<
      stomp_map->RegionArea(region_iter) << " sq. degrees.\n";

  // Now we include a whole series of RegionBounds that should divide up the
  // Map area with nothing left over.  The RegionBounds are also set up so
  // that they have unequal areas so that we test the ability of the code to
  // handle that case.
  regionVec.clear();

  Stomp::LatLonBound* upper_left =
    new Stomp::LatLonBound(0.0, 5.0, -5.0, 0.0,
			   Stomp::AngularCoordinate::Equatorial);
  Stomp::RegionBound upper_left_region(upper_left);
  regionVec.push_back(upper_left_region);

  Stomp::LatLonBound* lower_left =
    new Stomp::LatLonBound(-5.0, 0.0, -5.0, 0.0,
			   Stomp::AngularCoordinate::Equatorial);
  Stomp::RegionBound lower_left_region(lower_left);
  regionVec.push_back(lower_left_region);

  Stomp::LatLonBound* upper_right =
    new Stomp::LatLonBound(0.0, 5.0, 0.0, 5.0,
			   Stomp::AngularCoordinate::Equatorial);
  Stomp::RegionBound upper_right_region(upper_right);
  regionVec.push_back(upper_right_region);

  Stomp::LatLonBound* lower_right0 =
    new Stomp::LatLonBound(-5.0, -4.0, 0.0, 5.0,
			   Stomp::AngularCoordinate::Equatorial);
  Stomp::RegionBound lower_right0_region(lower_right0);
  regionVec.push_back(lower_right0_region);

  Stomp::LatLonBound* lower_right1 =
    new Stomp::LatLonBound(-4.0, -2.5, 0.0, 5.0,
			   Stomp::AngularCoordinate::Equatorial);
  Stomp::RegionBound lower_right1_region(lower_right1);
  regionVec.push_back(lower_right1_region);

  Stomp::LatLonBound* lower_right2 =
    new Stomp::LatLonBound(-2.5, 0.0, 0.0, 5.0,
			   Stomp::AngularCoordinate::Equatorial);
  Stomp::RegionBound lower_right2_region(lower_right2);
  regionVec.push_back(lower_right2_region);

  std::cout << "\nRe-Regionating the map with multiple RegionBounds...\n";
  n_regions_final = stomp_map->InitializeRegions(regionVec, n_regions);
  std::cout << "\tRegionated into " << n_regions_final << " (" <<
    n_regions << ") regions at " << stomp_map->RegionResolution() << "\n";
  for (uint16_t region_iter=0;region_iter<n_regions_final;region_iter++)
    std::cout << "\t\t" << region_iter << ": " <<
      stomp_map->RegionArea(region_iter) << " sq. degrees.\n";
}

void MapSoftenTests() {
  // Ok, now we test our routines for softening the maximum resolution of the
  // Map and cut the map based on the Weight.
  std::cout << "\n";
  std::cout << "************************\n";
  std::cout << "*** Map Soften Tests ***\n";
  std::cout << "************************\n\n";
  double theta = 3.0;
  uint32_t resolution = 256;
  uint32_t soft_resolution = resolution/8;
  Stomp::AngularCoordinate ang(60.0, 0.0, Stomp::AngularCoordinate::Survey);
  Stomp::Pixel tmp_pix(ang, resolution);
  Stomp::PixelVector annulus_pix;
  tmp_pix.WithinRadius(theta, annulus_pix);
  Stomp::Map* stomp_map = new Stomp::Map(annulus_pix);

  // Before softening, we check some global properties of the starting Map.
  std::cout << "Starting Area: " << stomp_map->Area() << " sq. deg.\n";
  std::cout << "\tMax resolution: " << stomp_map->MaxResolution() <<
    ", Min resolution: " << stomp_map->MinResolution() << "\n";
  std::cout << "\tMax weight: " << stomp_map->MaxWeight() <<
    ", Min weight: " << stomp_map->MinWeight() << "\n";

  // Now we create a softened version of the starting map with average weights.
  Stomp::Map soft_map;
  std::cout << "\nSoftening to " << soft_resolution <<
    " with average weights...\n";
  stomp_map->Soften(soft_map, soft_resolution, true);
  std::cout << "\tArea: " << soft_map.Area() << " sq. deg.\n";
  std::cout << "\tMax resolution: " << soft_map.MaxResolution() <<
    ", Min resolution: " << soft_map.MinResolution() << "\n";
  std::cout << "\tMax weight: " << soft_map.MaxWeight() <<
    ", Min weight: " << soft_map.MinWeight() << "\n";

  // Now we create a softened version of the starting map with unity weights.
  std::cout << "\nSoftening to " << soft_resolution <<
    " with unmasked fractions...\n";
  stomp_map->Soften(soft_map, soft_resolution, false);
  std::cout << "\tArea: " << soft_map.Area() << " sq. deg.\n";
  std::cout << "\tMax resolution: " << soft_map.MaxResolution() <<
    ", Min resolution: " << soft_map.MinResolution() << "\n";
  std::cout << "\tMax weight: " << soft_map.MaxWeight() <<
    ", Min weight: " << soft_map.MinWeight() << "\n";

  // Now let's soften the current map with unity weights and make some cuts
  // on weight.
  double weight_limit = 0.51;
  stomp_map->Soften(soft_resolution, false);
  std::cout << "\nSplitting map at Weight() = " << weight_limit << "...\n";
  stomp_map->SetMaximumWeight(weight_limit);
  std::cout << "\nSoftening weight <= " << weight_limit << "...\n";
  std::cout << "\tArea: " << stomp_map->Area() << " sq. deg.\n";
  std::cout << "\tMax resolution: " << stomp_map->MaxResolution() <<
    ", Min resolution: " << stomp_map->MinResolution() << "\n";
  std::cout << "\tMax weight: " << stomp_map->MaxWeight() <<
    ", Min weight: " << stomp_map->MinWeight() << "\n";

  soft_map.SetMinimumWeight(weight_limit);
  std::cout << "\nSoftening weight >= " << weight_limit << "...\n";
  std::cout << "\tArea: " << soft_map.Area() << " sq. deg.\n";
  std::cout << "\tMax resolution: " << soft_map.MaxResolution() <<
    ", Min resolution: " << soft_map.MinResolution() << "\n";
  std::cout << "\tMax weight: " << soft_map.MaxWeight() <<
    ", Min weight: " << soft_map.MinWeight() << "\n";
}

// Define our command line flags here so we can use these flags later.
DEFINE_bool(all_map_tests, false, "Run all class unit tests.");
DEFINE_bool(map_basic_tests, false, "Run Map basic tests");
DEFINE_bool(map_write_tests, false, "Run Map write tests");
DEFINE_bool(map_read_tests, false, "Run Map read tests");
DEFINE_bool(map_pixelization_tests, false, "Run Map pixelization tests");
DEFINE_bool(map_cover_tests, false, "Run Map cover tests");
DEFINE_bool(map_iterator_tests, false, "Run Map iterator tests");
DEFINE_bool(map_location_tests, false, "Run Map location tests");
DEFINE_bool(map_unmasked_fraction_tests, false,
            "Run Map unmasked fraction tests");
DEFINE_bool(map_contains_tests, false, "Run Map Contains tests");
DEFINE_bool(map_random_points_tests, false, "Run Map random points tests");
DEFINE_bool(map_multimap_tests, false, "Run Map multi-map tests");
DEFINE_bool(map_region_tests, false, "Run Map region tests");
DEFINE_bool(map_region_bound_tests, false, "Run Map RegionBound tests");
DEFINE_bool(map_soften_tests, false, "Run Map soften tests");

void MapUnitTests(bool run_all_tests) {
  void MapBasicTests();
  void MapWriteTests();
  void MapReadTests();
  void MapPixelizationTests();
  void MapCoverTests();
  void MapIteratorTests();
  void MapLocationTests();
  void MapUnmaskedFractionTests();
  void MapContainsTests();
  void MapRandomPointsTests();
  void MapMultiMapTests();
  void MapRegionTests();
  void MapRegionBoundTests();
  void MapSoftenTests();

  if (run_all_tests) FLAGS_all_map_tests = true;

  // Check the basic routines for generating Stomp::Map instances from a list
  // of contiguous pixels.
  if (FLAGS_all_map_tests || FLAGS_map_basic_tests) MapBasicTests();

  // Check the routines for writing Stomp::Map's to file.
  if (FLAGS_all_map_tests || FLAGS_map_write_tests) MapWriteTests();

  // Check the routines for reading a Stomp::Map from a simple ASCII file.
  if (FLAGS_all_map_tests || FLAGS_map_read_tests) MapReadTests();

  // Check the routines for generating a Stomp::Map from a GeometricBound
  if (FLAGS_all_map_tests || FLAGS_map_pixelization_tests)
    MapPixelizationTests();

  // Check the routines for finding the Stomp::Map Coverage and Covering.
  if (FLAGS_all_map_tests || FLAGS_map_cover_tests) MapCoverTests();

  // Check the routines for iterating through a Stomp::Map's pixels using the
  // Begin, End and Iterate methods.
  if (FLAGS_all_map_tests || FLAGS_map_iterator_tests) MapIteratorTests();

  // Check the Stomp::Map methods for checking locations against the area
  // covered by a Stomp::Map.
  if (FLAGS_all_map_tests || FLAGS_map_location_tests) MapLocationTests();

  // Check the routines for finding the area of a pixel covered by a Stomp::Map.
  if (FLAGS_all_map_tests || FLAGS_map_unmasked_fraction_tests)
    MapUnmaskedFractionTests();

  // Check the routines for testing whether one Stomp::Map is contained in
  // another Stomp::Map.
  if (FLAGS_all_map_tests || FLAGS_map_contains_tests) MapContainsTests();

  // Check the routine for generating random points within a Stomp::Map's area.
  if (FLAGS_all_map_tests || FLAGS_map_random_points_tests)
    MapRandomPointsTests();

  // Check the different ways of combining Stomp::Map instances.
  if (FLAGS_all_map_tests || FLAGS_map_multimap_tests) MapMultiMapTests();

  // Check the routines for breaking a Stomp::Map into sub-regions.
  if (FLAGS_all_map_tests || FLAGS_map_region_tests) MapRegionTests();

  // Check the routines for breaking a Stomp::Map into sub-regions with
  // specified RegionBounds.
  if (FLAGS_all_map_tests || FLAGS_map_region_bound_tests)
    MapRegionBoundTests();

  // Check the routines for softening the maximum resolution of the
  // Map and cutting the map based on the Weight.
  if (FLAGS_all_map_tests || FLAGS_map_soften_tests) MapSoftenTests();
}
