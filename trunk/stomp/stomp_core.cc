// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This file contains the basic constants necessary to make the library
// function and should be included in any program using the library.


#include <stomp/stomp_core.h>
#include <math.h>

namespace Stomp {

// Define our Stomp class constants.  See stomp_util.h for their meanings.
const double Pi = 2.0*asin(1.0);
const double DegToRad = Pi/180.0;
const double RadToDeg = 180.0/Pi;
const double StradToDeg = 180.0*180.0/(Pi*Pi);
const uint32_t Nx0 = 36;
const uint32_t Ny0 = 13;
const double EtaOffSet = 91.25;
const double SurveyCenterRA = 185.0;
const double SurveyCenterDEC = 32.5;
const double Node = DegToRad*(SurveyCenterRA-90.0);
const double EtaPole = DegToRad*SurveyCenterDEC;
const uint8_t HPixLevel = 2;  // 4 in the old parlance.
const uint8_t MaxPixelLevel = 15;  // 2^15 = 32768
const uint32_t HPixResolution = 1 << HPixLevel;  // 4
const uint32_t MaxPixelResolution = 1 << MaxPixelLevel;  // 32768
const uint8_t ResolutionLevels = MaxPixelLevel - HPixLevel + 1;
const double HPixArea =
      4.0*Pi*StradToDeg/(HPixResolution*HPixResolution*Nx0*Ny0);
const uint32_t MaxPixnum = Nx0*Ny0*2048*2048;
const uint32_t MaxSuperpixnum = Nx0*Ny0*HPixResolution*HPixResolution;

bool DoubleLT(double a, double b) {
  return (a < b - 1.0e-10 ? true : false);
}

bool DoubleLE(double a, double b) {
  return (a <= b + 1.0e-10 ? true : false);
}

bool DoubleGT(double a, double b) {
  return (a > b + 1.0e-10 ? true : false);
}

bool DoubleGE(double a, double b) {
  return (a >= b - 1.0e-10 ? true : false);
}

bool DoubleEQ(double a, double b) {
  return (DoubleLE(a, b) && DoubleGE(a, b) ? true : false);
}

uint8_t MostSignificantBit(uint32_t input_int) {
  uint8_t ln_int = 0;
  while (input_int >>= 1) ln_int++;
  return ln_int;
}

}  // end namespace Stomp
