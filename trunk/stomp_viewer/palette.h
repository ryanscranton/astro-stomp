// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the Palette class.  It handles the translation of
// Stomp::Map and Stomp::WAngularVector weights into colors for display.

#ifndef PALETTE_H
#define PALETTE_H

#include <QColor>
#include <stomp.h>

class Palette;

class Palette {
  // This class handles the translation between map or point weight and color.
  // The PaletteType enum determines the color scheme.  In general, the color
  // returned by the weight is a linear interpolation along the palette between
  // the weight minimum and maximum.  A logrithmic scaling is available if that
  // option is selected.

 public:
  enum PaletteType {
    BlueTemperature,
    GreenTemperature,
    RedTemperature,
    GrayScale,
    InverseGrayScale,
    Rainbow,
    InverseRainbow
  };
  Palette();
  Palette(PaletteType palette, double weight_min, double weight_max,
	  bool log_weights);
  ~Palette();
  void initialize(PaletteType palette);
  void initializeBlueTemperature();
  void initializeGreenTemperature();
  void initializeRedTemperature();
  void initializeGrayScale();
  void initializeInverseGrayScale();
  void initializeRainbow();
  void initializeInverseRainbow();
  void setWeightRange(double weight_min, double weight_max);
  void useLogWeights(bool use_log_weights);
  double weightMin();
  double weightMax();
  bool logWeight();
  QRgb rgb(double weight);
  QColor color(double weight);
  int colorIndex(double weight);
  std::string currentPalette();
  PaletteType currentPaletteType();

 private:
  PaletteType palette_type_;
  std::vector<QRgb> rgb_;
  double weight_min_, weight_max_;
  bool log_weight_;
};

#endif
