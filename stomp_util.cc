#include <math.h>
#include <iostream>
#include "stomp_core.h"
#include "stomp_util.h"

namespace Stomp {

// Define our default WMAP5 flat, LCDM cosmology.
double Cosmology::omega_m = 0.2736;
double Cosmology::h = 0.705;
const double Cosmology::AA_ = 1.718;
const double Cosmology::BB_ = 0.315;
double Cosmology::a_ = Cosmology::AA_*Cosmology::omega_m;
double Cosmology::b_ = Cosmology::BB_*sqrt(Cosmology::omega_m);

double Cosmology::OmegaM() {
  return omega_m;
}

double Cosmology::HubbleConstant() {
  return h*100.0;
}

double Cosmology::HubbleDistance() {
  return 3000.0/h;
}

double Cosmology::OmegaL() {
  return 1.0 - omega_m;
}

void Cosmology::SetOmegaM(double new_omega_m) {
  omega_m = new_omega_m;
  a_ = AA_*omega_m;
  b_ = BB_*sqrt(omega_m);
}

void Cosmology::SetHubbleConstant(double hubble) {
    h = hubble/100.0;
}

void Cosmology::SetOmegaL(double omega_lambda) {
  omega_m = 1.0 - omega_lambda;
  a_ = AA_*omega_m;
  b_ = BB_*sqrt(omega_m);
}

double Cosmology::ComovingDistance(double z) {
  // In Mpc/h.
  return HubbleDistance()*z/sqrt(1.0 + a_*z + b_*z*z);
}

double Cosmology::AngularDiameterDistance(double z) {
  // In Mpc/h.
  return ComovingDistance(z)/(1.0+z);
}

double Cosmology::LuminosityDistance(double z) {
  // In Mpc/h.
  return ComovingDistance(z)*(1.0+z);
}

double Cosmology::ProjectedDistance(double z, double theta) {
  // where theta is assumed to be in degrees and the return value in Mpc/h.
  return theta*DegToRad*AngularDiameterDistance(z);
}

double Cosmology::ProjectedAngle(double z, double radius) {
  // where radius is in Mpc/h and the return angle is in degrees.
  return RadToDeg*radius/AngularDiameterDistance(z);
}

StompWatch::StompWatch() {
  gettimeofday(&start, NULL);
  stop = start;
}

void StompWatch::StartTimer() {
  gettimeofday(&start, NULL);
}

void StompWatch::StopTimer() {
  gettimeofday(&stop, NULL);
}

double StompWatch::ElapsedTime() {
  timeval interval;
  timersub(&stop, &start, &interval);
  // Need to convert tv_usec from microseconds to seconds.
  return interval.tv_sec + interval.tv_usec/1000000.0;
}

} // end namespace Stomp
