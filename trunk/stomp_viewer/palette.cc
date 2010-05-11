#include <QtGui>
#include "palette.h"

Palette::Palette() {
  initialize(BlueTemperature);
  weight_min_ = 0.0;
  weight_max_ = 1.0;
  log_weight_ = false;
}

Palette::Palette(PaletteType palette_type, double weight_min,
		 double weight_max, bool log_weight) {
  initialize(palette_type);
  setWeightRange(weight_min, weight_max);
  useLogWeights(log_weight);
}

Palette::~Palette() {
  rgb_.clear();
}

void Palette::initialize(PaletteType palette_type) {
  if (!rgb_.empty()) rgb_.clear();
  rgb_.reserve(256);

  palette_type_ = palette_type;

  switch (palette_type_) {
  case BlueTemperature:
    initializeBlueTemperature();
    break;
  case GreenTemperature:
    initializeGreenTemperature();
    break;
  case RedTemperature:
    initializeRedTemperature();
    break;
  case GrayScale:
    initializeGrayScale();
    break;
  case InverseGrayScale:
    initializeInverseGrayScale();
    break;
  case Rainbow:
    initializeRainbow();
    break;
  case InverseRainbow:
    initializeInverseRainbow();
    break;
  }
}

void Palette::initializeBlueTemperature() {
  for (int i=0;i<256;i++) rgb_[i] = QColor(0, 0, i).rgba();
}

void Palette::initializeGreenTemperature() {
  for (int i=0;i<256;i++) rgb_[i] = QColor(0, i, 0).rgba();
}

void Palette::initializeRedTemperature() {
  for (int i=0;i<256;i++) rgb_[i] = QColor(i, 0, 0).rgba();
}

void Palette::initializeGrayScale() {
  for (int i=0;i<256;i++) rgb_[i] = QColor(i, i, i).rgba();
}

void Palette::initializeInverseGrayScale() {
  for (int i=0;i<256;i++)
    rgb_[i] = QColor(255 - i, 255 - i, 255 - i).rgba();
}

void Palette::initializeRainbow() {
  double step = 300.0/255.0;
  for (int i=0;i<256;i++) {
    QColor qcolor;
    qcolor.setHsv(static_cast<int>(step*i), 255, 255);
    rgb_[i] = qcolor.rgba();
  }
}

void Palette::initializeInverseRainbow() {
  double step = 300.0/255.0;
  for (int i=0;i<256;i++) {
    QColor qcolor;
    qcolor.setHsv(static_cast<int>(300.0 - step*i), 255, 255);
    rgb_[i] = qcolor.rgba();
  }
}

void Palette::setWeightRange(double weight_min, double weight_max) {
  weight_min_ = weight_min;
  weight_max_ = weight_max;
}

void Palette::useLogWeights(bool use_log_weights) {
  log_weight_ = use_log_weights;

  if (log_weight_) {
    if (Stomp::DoubleLT(weight_min_, 0.0)) weight_min_ = 1.0e-10;
    if (Stomp::DoubleLT(weight_max_, weight_min_))
      weight_max_ = weight_min_ + 1.0e-10;
  }
}

double Palette::weightMin() {
  return weight_min_;
}

double Palette::weightMax() {
  return weight_max_;
}

bool Palette::logWeight() {
  return log_weight_;
}

QRgb Palette::rgb(double weight) {
  if (Stomp::DoubleGE(weight, weight_max_)) weight = weight_max_;
  if (Stomp::DoubleLE(weight, weight_min_)) weight = weight_min_;

  double norm_weight = (weight - weight_min_)/(weight_max_ - weight_min_);
  if (log_weight_) {
    norm_weight = log10(weight) - log10(weight_min_);
    norm_weight /= log10(weight_max_) - log10(weight_min_);
  }

  int i = static_cast<int>(256.0*norm_weight);
  if (i < 0) i = 0;
  if (i > 255) i = 255;
  return rgb_[i];
}

QColor Palette::color(double weight) {
  return QColor(rgb(weight));
}

int Palette::colorIndex(double weight) {
  if (Stomp::DoubleGE(weight, weight_max_)) weight = weight_max_;
  if (Stomp::DoubleLE(weight, weight_min_)) weight = weight_min_;

  double norm_weight = (weight - weight_min_)/(weight_max_ - weight_min_);
  if (log_weight_) {
    norm_weight = log10(weight) - log10(weight_min_);
    norm_weight /= log10(weight_max_) - log10(weight_min_);
  }

  int i = static_cast<int>(256.0*norm_weight);
  if (i < 0) i = 0;
  if (i > 255) i = 255;
  return i;
}

std::string Palette::currentPalette() {
  std::string current_palette;

  switch (palette_type_) {
  case BlueTemperature:
    current_palette = "BlueTemperature";
    break;
  case GreenTemperature:
    current_palette = "GreenTemperature";
    break;
  case RedTemperature:
    current_palette = "RedTemperature";
    break;
  case GrayScale:
    current_palette = "GrayScale";
    break;
  case InverseGrayScale:
    current_palette = "InverseGrayScale";
    break;
  case Rainbow:
    current_palette = "Rainbow";
    break;
  case InverseRainbow:
    current_palette = "InverseRainbow";
    break;
  }
  return current_palette;
}

Palette::PaletteType Palette::currentPaletteType() {
  return palette_type_;
}

