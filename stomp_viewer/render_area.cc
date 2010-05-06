#include <QtGui>
#include "render_area.h"

Palette::Palette() {
  Initialize(BlueTemperature);
}

Palette::Palette(PaletteType palette_type) {
  Initialize(palette_type);
}

Palette::~Palette() {
  rgb_.clear();
}

void Palette::Initialize(PaletteType palette_type) {
  if (rgb_.empty()) rgb_.reserve(256);
  palette_type_ = palette_type;

  switch (palette_type_) {
  case BlueTemperature:
    InitializeBlueTemperature();
    break;
  case GreenTemperature:
    InitializeGreenTemperature();
    break;
  case RedTemperature:
    InitializeRedTemperature();
    break;
  case GrayScale:
    InitializeGrayScale();
    break;
  case InverseGrayScale:
    InitializeInverseGrayScale();
    break;
  case Rainbow:
    InitializeRainbow();
    break;
  case InverseRainbow:
    InitializeInverseRainbow();
    break;
  }
}

void Palette::InitializeBlueTemperature() {
  for (uint16_t i=0;i<256;i++) rgb_[i] = QColor(0, 0, i);
}

void Palette::InitializeGreenTemperature() {
  for (uint16_t i=0;i<256;i++) rgb_[i] = QColor(0, i, 0);
}

void Palette::InitializeRedTemperature() {
  for (uint16_t i=0;i<256;i++) rgb_[i] = QColor(i, 0, 0);
}

void Palette::InitializeGrayScale() {
  for (uint16_t i=0;i<256;i++) rgb_[i] = QColor(i, i, i);
}

void Palette::InitializeInverseGrayScale() {
  for (uint16_t i=0;i<256;i++) rgb_[i] = QColor(255 - i, 255 - i, 255 - i);
}

void Palette::InitializeRainbow() {
  double step = 300.0/255.0;
  for (uint16_t i=0;i<256;i++) {
    QColor qcolor;
    qcolor.setHsv(static_cast<int>(step*i), 255, 255);
    rgb_[i] = qcolor;
  }
}

void Palette::InitializeInverseRainbow() {
  double step = 300.0/255.0;
  for (uint16_t i=0;i<256;i++) {
    QColor qcolor;
    qcolor.setHsv(static_cast<int>(300.0 - step*i), 255, 255);
    rgb_[i] = qcolor;
  }
}

QColor Palette::Color(double weight) {
  int i = static_cast<int>(256.0*weight);
  if (i < 0) i = 0;
  if (i > 255) i = 255;
  return rgb_[i];
}

std::string Palette::CurrentPalette() {
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
Palette::PaletteType Palette::CurrentPaletteType() {
  return palette_type_;
}

RenderArea::RenderArea(QWidget *parent) : QWidget(parent) {
  antialiased_ = false;
  auto_update_ = false;

  setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

  palette_ = Palette();
  points_palette_ = Palette();

  Stomp::PixelVector pix;
  double dweight = 1.0/(4*Stomp::MaxSuperpixnum);
  for (uint32_t k=0;k<4*Stomp::MaxSuperpixnum;k++)
    pix.push_back(Stomp::Pixel(8, k, k*dweight));
  Stomp::Map* base_map = new Stomp::Map(pix);

  for (uint8_t level=Stomp::HPixLevel;level<=Stomp::MaxPixelLevel;level++)
    stomp_map_[level] = base_map;
  good_map_ = true;
  good_points_ = false;

  buffer_top_ = buffer_bottom_ = buffer_left_ = buffer_right_ = 0;
  fill_ = false;

  sphere_ = Stomp::AngularCoordinate::Equatorial;
  setFullSky(false);
  useAitoffProjection(false);

  weight_min_ = 0.0;
  weight_max_ = 1.0;
  points_weight_min_ = 0.0;
  points_weight_max_ = 0.0;

  max_resolution_ = Stomp::HPixResolution;

  setBackgroundRole(QPalette::Base);
  setAutoFillBackground(true);

  background_color_ = Qt::white;
  pixmap_ = new QPixmap(width(), height());
  pixmap_->fill(background_color_);

  show_grid_ = false;
  show_coordinates_ = false;
  grid_pixmap_ = new QPixmap(width(), height());
  grid_pixmap_->fill(Qt::transparent);

  show_points_ = false;
  filter_points_ = false;
  points_pixmap_ = new QPixmap(width(), height());
  points_pixmap_->fill(Qt::transparent);

  n_tick_divisions_ = 21;
  tick_divisions_ = new double[n_tick_divisions_];
  tick_divisions_[0] = 0.00001;
  tick_divisions_[1] = 0.00002;
  tick_divisions_[2] = 0.00005;
  tick_divisions_[3] = 0.0001;
  tick_divisions_[4] = 0.0002;
  tick_divisions_[5] = 0.0005;
  tick_divisions_[6] = 0.001;
  tick_divisions_[7] = 0.002;
  tick_divisions_[8] = 0.005;
  tick_divisions_[9] = 0.01;
  tick_divisions_[10] = 0.02;
  tick_divisions_[11] = 0.05;
  tick_divisions_[12] = 0.1;
  tick_divisions_[13] = 0.2;
  tick_divisions_[14] = 0.5;
  tick_divisions_[15] = 1.0;
  tick_divisions_[16] = 2.0;
  tick_divisions_[17] = 5.0;
  tick_divisions_[18] = 10.0;
  tick_divisions_[19] = 15.0;
  tick_divisions_[20] = 30.0;

  major_tick_length_ = 0;
  minor_tick_length_ = 0;
  tick_precision_ = 4;
  minor_ticks_per_major_ = 4;
}

QSize RenderArea::minimumSizeHint() const {
  return QSize(200, 100);
}

QSize RenderArea::sizeHint() const {
  return QSize(1024, 512);
}

double RenderArea::weightMin() {
  return weight_min_;
}

double RenderArea::weightMax() {
  return weight_max_;
}

double RenderArea::pointsWeightMin() {
  return points_weight_min_;
}

double RenderArea::pointsWeightMax() {
  return points_weight_max_;
}

double RenderArea::longitudeMin() {
  return lonmin_;
}

double RenderArea::longitudeMax() {
  return lonmax_;
}

double RenderArea::latitudeMin() {
  return latmin_;
}

double RenderArea::latitudeMax() {
  return latmax_;
}

uint32_t RenderArea::maxResolution() {
  return max_resolution_;
}

bool RenderArea::autoUpdating() {
  return auto_update_;
}

bool RenderArea::fullSky() {
  return full_sky_;
}

bool RenderArea::displayingCoordinates() {
  return show_coordinates_;
}

bool RenderArea::displayingGrid() {
  return show_grid_;
}

bool RenderArea::displayingPoints() {
  return show_points_;
}

bool RenderArea::filteringPoints() {
  return filter_points_;
}

void RenderArea::updatePixmap(bool send_update_call) {
  delete pixmap_;

  pixmap_ = new QPixmap(size());
  pixmap_->fill(background_color_);

  paintMap(pixmap_);
  updatePoints(false);
  updateGrid(false);
  if (send_update_call) update();
}

void RenderArea::updatePoints(bool send_update_call) {
  delete points_pixmap_;

  points_pixmap_ = new QPixmap(size());
  points_pixmap_->fill(Qt::transparent);

  paintPoints(points_pixmap_);
  if (send_update_call) update();
}

void RenderArea::updateGrid(bool send_update_call) {
  delete grid_pixmap_;

  grid_pixmap_ = new QPixmap(size());
  grid_pixmap_->fill(Qt::transparent);

  paintGrid(grid_pixmap_);
  if (send_update_call) update();
}

bool RenderArea::readNewMap(QString& input_file) {
  _clearMaps();

  std::string stl_input_file = std::string(input_file.toLatin1());
  Stomp::Map* base_map = new Stomp::Map(stl_input_file);

  good_map_ = !base_map->Empty();

  if (good_map_) {
    _findNewWeightBounds(base_map);
    _findNewImageBounds(base_map);
    _findNewMaxResolution(base_map);
    emit newMapParameters();


    // Now, we iterate over our possible resolution values for rendering
    // the Map.  If the resolution limit is higher than the maximum resolution
    // of the basic Map, then we put assign a copy of the pointer to
    // that resolution.
    for (uint8_t level=Stomp::MaxPixelLevel;level>=Stomp::HPixLevel;level--) {
      uint32_t resolution = 1 << level;
      if (resolution >= max_resolution_) {
	stomp_map_[level] = base_map;
      } else {
	// For the lower resolution versions of the map, we create softened
	// versions of the map.  We can do this progressively, so that the
	// computation is faster than resampling from the basic Map every time.
	Stomp::Map* soft_map = new Stomp::Map();

	stomp_map_[level+1]->Soften(*soft_map, resolution, true);
	stomp_map_[level] = soft_map;
      }
    }

    // When we render the map for the first time, we want to plot enough pixels
    // to give a good idea of what the map looks like, but avoid spending
    // several minutes rendering small pixels that won't actually show up.
    // To balance this out, we choose an intial rendering resolution where
    // the total number of pixels is a few thousand.
    _findRenderLevel(8000);

    updatePixmap();
  }

  return good_map_;
}

bool RenderArea::readNewPointsFile(QString& input_file,
				   Stomp::AngularCoordinate::Sphere sphere,
				   bool weighted_points) {
  ang_.clear();

  std::string stl_input_file = std::string(input_file.toLatin1());

  std::ifstream input_file_ptr(stl_input_file.c_str());

  bool good_points_ = false;

  if (input_file_ptr) {
    good_points_ = true;

    double theta, phi, weight = 1.0;

    while (!input_file_ptr.eof()) {
      if (weighted_points) {
	input_file_ptr >> theta >> phi >> weight;
      } else {
	input_file_ptr >> theta >> phi;
      }

      if (!input_file_ptr.eof()) {
	Stomp::WeightedAngularCoordinate tmp_ang(theta, phi, 1.0, sphere);
	ang_.push_back(tmp_ang);
      }
    }
    input_file_ptr.close();

    _findNewWeightBounds(ang_);
    if (!good_map_) _findNewImageBounds(ang_);
    emit newMapParameters();

    updatePoints();
  }

  return good_points_;
}

void RenderArea::clearMap() {
  _clearMaps();
  updatePixmap();
}

void RenderArea::clearPoints() {
  ang_.clear();
  updatePoints();
}

bool RenderArea::writeToPng(QString& output_file) {
  if (size() != pixmap_->size()) updatePixmap();

  bool painted_map = pixmap_->save(output_file);

  return painted_map;
}

void RenderArea::setAntialiased(bool antialiased) {
  this->antialiased_ = antialiased;
  if (auto_update_) updatePixmap();
}

void RenderArea::setAutoUpdate(bool auto_update) {
  this->auto_update_ = auto_update;
}

void RenderArea::setWeightRange(double weight_min, double weight_max) {
  this->weight_min_ = weight_min;
  this->weight_max_ = weight_max;
  if (auto_update_) updatePixmap();
}

void RenderArea::setPointsWeightRange(double weight_min, double weight_max) {
  this->points_weight_min_ = weight_min;
  this->points_weight_max_ = weight_max;
  if (auto_update_) updatePixmap();
}

void RenderArea::setPalette(Palette::PaletteType palette_type) {
  this->palette_.Initialize(palette_type);
  if (auto_update_) updatePixmap();
}

void RenderArea::setPointsPalette(Palette::PaletteType palette_type) {
  this->points_palette_.Initialize(palette_type);
  if (auto_update_) updatePixmap();
}

void RenderArea::setVerticalBuffer(uint16_t buffer_pixels) {
  buffer_top_ = buffer_pixels;
  buffer_bottom_ = buffer_pixels;
  if (auto_update_) updatePixmap();
}

void RenderArea::setHorizontalBuffer(uint16_t buffer_pixels) {
  buffer_left_ = buffer_pixels;
  buffer_right_ = buffer_pixels;
  if (auto_update_) updatePixmap();
}

void RenderArea::setTopBuffer(uint16_t buffer_pixels) {
  buffer_top_ = buffer_pixels;
  if (auto_update_) updatePixmap();
}

void RenderArea::setBottomBuffer(uint16_t buffer_pixels) {
  buffer_bottom_ = buffer_pixels;
  if (auto_update_) updatePixmap();
}

void RenderArea::setLeftBuffer(uint16_t buffer_pixels) {
  buffer_left_ = buffer_pixels;
  if (auto_update_) updatePixmap();
}

void RenderArea::setRightBuffer(uint16_t buffer_pixels) {
  buffer_right_ = buffer_pixels;
  if (auto_update_) updatePixmap();
}

void RenderArea::setImageBounds(double lonmin, double lonmax,
				double latmin, double latmax) {
  latmin_ = latmin;
  latmax_ = latmax;
  lonmin_ = lonmin;
  lonmax_ = lonmax;
  _checkBounds();
  if (auto_update_) updatePixmap();
}

void RenderArea::setImageBounds(double lonmin, double lonmax,
				double latmin, double latmax,
				Stomp::AngularCoordinate::Sphere sphere) {
  latmin_ = latmin;
  latmax_ = latmax;
  lonmin_ = lonmin;
  lonmax_ = lonmax;
  sphere_ = sphere;
  _checkBounds();
  if (auto_update_) updatePixmap();
}

void RenderArea::setSurveyBounds(double lammin, double lammax,
				 double etamin, double etamax) {
  latmin_ = lammin;
  latmax_ = lammax;
  lonmin_ = etamin;
  lonmax_ = etamax;
  sphere_ = Stomp::AngularCoordinate::Survey;
  _checkBounds();
  if (auto_update_) updatePixmap();
}

void RenderArea::setEquatorialBounds(double ramin, double ramax,
				     double decmin, double decmax) {
  latmin_ = decmin;
  latmax_ = decmax;
  lonmin_ = ramin;
  lonmax_ = ramax;
  sphere_ = Stomp::AngularCoordinate::Equatorial;
  _checkBounds();
  if (auto_update_) updatePixmap();
}

void RenderArea::setGalacticBounds(double lonmin, double lonmax,
				   double latmin, double latmax) {
  latmin_ = latmin;
  latmax_ = latmax;
  lonmin_ = lonmin;
  lonmax_ = lonmax;
  sphere_ = Stomp::AngularCoordinate::Galactic;
  _checkBounds();
  if (auto_update_) updatePixmap();
}

void RenderArea::setFullSky(bool full_sky) {
  if (full_sky) {
    latmin_ = -90.0;
    latmax_ = 90.0;
    if (sphere_ == Stomp::AngularCoordinate::Survey) {
      lonmin_ = -180.0;
      lonmax_ = 180.0;
    } else {
      lonmin_ = 0.0;
      lonmax_ = 360.0;
    }
    full_sky_ = true;
    continuous_longitude_ = true;
    if (auto_update_) updatePixmap();
  } else {
    full_sky_ = false;
  }
}

void RenderArea::setCoordinateSystem(Stomp::AngularCoordinate::Sphere sphere) {
  sphere_ = sphere;
  if (full_sky_) {
    setFullSky(true);
  } else {
    _findNewImageBounds(stomp_map_[Stomp::HPixResolution]);
  }
  if (auto_update_) updatePixmap();
}

void RenderArea::useSurveyCoordinates() {
  sphere_ = Stomp::AngularCoordinate::Survey;
  setFullSky(full_sky_);
  _checkBounds();
  if (auto_update_) updatePixmap();
}

void RenderArea::useEquatorialCoordinates() {
  sphere_ = Stomp::AngularCoordinate::Equatorial;
  setFullSky(full_sky_);
  _checkBounds();
  if (auto_update_) updatePixmap();
}

void RenderArea::useGalacticCoordinates() {
  sphere_ = Stomp::AngularCoordinate::Galactic;
  setFullSky(full_sky_);
  _checkBounds();
  if (auto_update_) updatePixmap();
}

void RenderArea::zoomIn() {
  zoomIn(_lonCenter(), _latCenter());
}

void RenderArea::zoomOut() {
  zoomOut(_lonCenter(), _latCenter());
}

void RenderArea::zoomIn(double zoom_factor) {
  zoomIn(_lonCenter(), _latCenter(), zoom_factor);
}

void RenderArea::zoomOut(double zoom_factor) {
  zoomOut(_lonCenter(), _latCenter(), zoom_factor);
}

void RenderArea::fillPolygons(bool fill_polygons) {
  fill_ = fill_polygons;
  if (auto_update_) updatePixmap();
}

void RenderArea::useAitoffProjection(bool aitoff_projection) {
  aitoff_projection_ = aitoff_projection;
  if (auto_update_) updatePixmap();
}

void RenderArea::setDisplayCoordinates(bool display_coordinates) {
  show_coordinates_ = display_coordinates;
  if (auto_update_) updateGrid();
}

void RenderArea::setDisplayGrid(bool display_grid) {
  show_grid_ = display_grid;
  if (auto_update_) updateGrid();
}

void RenderArea::setDisplayPoints(bool display_points) {
  show_points_ = display_points;
  if (auto_update_) updatePoints();
}

void RenderArea::setFilterPoints(bool filter_points) {
  filter_points_ = filter_points;
  if (auto_update_) updatePoints();
}

void RenderArea::setMaxResolution(uint32_t resolution) {
  max_resolution_ = resolution;
  if (auto_update_) updatePixmap();
}

void RenderArea::setRenderResolution(uint32_t resolution) {
  render_level_ = Stomp::MostSignificantBit(resolution);
  if (auto_update_) updatePixmap();
}

void RenderArea::paintEvent(QPaintEvent * /* event */) {
  if (size() != pixmap_->size()) updatePixmap(false);

  QPainter painter(this);
  painter.drawPixmap(0, 0, *pixmap_);
  if (show_points_) painter.drawPixmap(0, 0, *points_pixmap_);
  if (show_grid_ || show_coordinates_) painter.drawPixmap(0, 0, *grid_pixmap_);
}

void RenderArea::mouseDoubleClickEvent(QMouseEvent *event) {
  QPointF new_center = event->posF();

  if (!aitoff_projection_)
    zoomIn(_xToLon(new_center.x()), _yToLat(new_center.y()));
}

void RenderArea::zoomIn(double new_longitude_center,
			double new_latitude_center,
			double zoom_factor) {
  // After zooming in, the longitude range (2*lon_span), should decrease by
  // 1/zoom_factor.
  double lon_span = 0.5*_lonRange()/zoom_factor;

  lonmin_ = new_longitude_center - lon_span;
  lonmax_ = new_longitude_center + lon_span;

  double lat_span = 0.5*_latRange()/zoom_factor;

  latmin_ = new_latitude_center - lat_span;
  latmax_ = new_latitude_center + lat_span;

  _checkBounds();

  updatePixmap();

  emit newMapParameters();
}

void RenderArea::zoomOut(double new_longitude_center,
			 double new_latitude_center,
			 double zoom_factor) {
  // After zooming out, the longitude range (2*lon_span), should increase by
  // zoom_factor.
  double lon_span = 0.5*_lonRange()*zoom_factor;
  double lat_span = 0.5*_latRange()*zoom_factor;

  if (Stomp::DoubleLE(lon_span, 180.0) &&
      Stomp::DoubleLE(lat_span, 90.0)) {
    // Provided that zooming out won't max out either of our coordinates, we
    // can safely increase each range.
    lonmin_ = new_longitude_center - lon_span;
    lonmax_ = new_longitude_center + lon_span;

    latmin_ = new_latitude_center - lat_span;
    latmax_ = new_latitude_center + lat_span;
  } else {
    if (Stomp::DoubleGE(lon_span, 180.0) &&
	Stomp::DoubleGE(lat_span, 90.0)) {
      // If increasing both ranges means we'd encompass at least the full sphere
      // then go to that option.
      setFullSky(true);
      setFullSky(false);
    } else {
      // Otherwise, we increase if we can or set to the maximum allowed range
      // if we can't.
      if (Stomp::DoubleLE(lon_span, 180.0)) {
	lonmin_ = new_longitude_center - lon_span;
	lonmax_ = new_longitude_center + lon_span;
      } else {
	if (sphere_ == Stomp::AngularCoordinate::Survey) {
	  lonmin_ = -180.0;
	  lonmax_ = 180.0;
	} else {
	  lonmin_ = 0.0;
	  lonmax_ = 360.0;
	}
      }
      if (Stomp::DoubleLE(lat_span, 90.0)) {
	latmin_ = new_latitude_center - lat_span;
	latmax_ = new_latitude_center + lat_span;
      } else {
	latmin_ = -90.0;
	latmax_ = 90.0;
      }
    }
  }

  _checkBounds();

  updatePixmap();

  emit newMapParameters();
}

bool RenderArea::paintMap(QPaintDevice *device) {
  QPainter painter;

  bool accessed_device = painter.begin(device);

  if (accessed_device) {
    if (antialiased_) {
      painter.setRenderHint(QPainter::Antialiasing, true);
      painter.translate(+0.5, +0.5);
    }

    if (!stomp_map_.empty()) {
      double render_pixel_area = _renderPixelArea();
      _findRenderLevel();
      Stomp::Map* render_map = stomp_map_[render_level_];

      for (Stomp::MapIterator iter=render_map->Begin();
	   iter!=render_map->End();render_map->Iterate(&iter)) {
	Stomp::PixelIterator pix = iter.second;
	if (pix->Resolution() <= max_resolution_ &&
	    (full_sky_ || pix->IntersectsBounds(lonmin_, lonmax_,
						latmin_, latmax_, sphere_))) {
	  QBrush pixel_brush =
	    QBrush(palette_.Color(_normalizeWeight(pix->Weight())));
	  painter.setBrush(pixel_brush);
	  painter.setPen(QPen(pixel_brush, 0));

	  if (pix->Area() < render_pixel_area) {
	    // If the pixel is so small that displaying it wouldn't cover at
	    // least a single pixel, we just use a single QPointF.
	    QPointF point;
	    _pixelToPoint(pix->PixelX(), pix->PixelY(), pix->Resolution(),
			  point);
	    painter.drawPoint(point);
	  } else {
	    if (pix->ContinuousBounds(sphere_)) {
	      QPolygonF polygon;
	      _pixelToPolygon(pix->PixelX(), pix->PixelY(),
			      pix->Resolution(), polygon);
	      if (fill_) {
		painter.drawConvexPolygon(polygon);
	      } else {
		painter.drawPolyline(polygon);
	      }
	    } else {
	      QPolygonF left_polygon;
	      QPolygonF right_polygon;
	      _splitPixelToPolygons(pix->PixelX(), pix->PixelY(),
				    pix->Resolution(), left_polygon,
				    right_polygon);
	      if (fill_) {
		painter.drawConvexPolygon(left_polygon);
		painter.drawConvexPolygon(right_polygon);
	      } else {
		painter.drawPolyline(left_polygon);
		painter.drawPolyline(right_polygon);
	      }
	    }
	  }
	}
      }
    }
    accessed_device = painter.end();
  }

  return accessed_device;
}

bool RenderArea::paintPoints(QPaintDevice *device) {
  QPainter painter;

  bool accessed_device = painter.begin(device);

  if (accessed_device) {
    if (antialiased_) {
      painter.setRenderHint(QPainter::Antialiasing, true);
      painter.translate(+0.5, +0.5);
    }

    if (!ang_.empty()) {
      for (uint32_t i=0;i<ang_.size();i++) {
	QPointF point;
	if (_angToPoint(ang_[i], point)) {
	  QBrush pixel_brush =
	    QBrush(points_palette_.
		   Color(_normalizePointsWeight(ang_[i].Weight())));
	  painter.setBrush(pixel_brush);
	  painter.setPen(QPen(pixel_brush, 0));
	  painter.drawPoint(point);
	}
      }
    }
    accessed_device = painter.end();
  }

  return accessed_device;
}

bool RenderArea::paintGrid(QPaintDevice *device) {
  QPainter painter;

  bool accessed_device = painter.begin(device);

  if (accessed_device) {
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.translate(+0.5, +0.5);

    painter.setPen(Qt::darkGray);

    QString lon_label;
    QString lat_label;

    switch (sphere_) {
    case Stomp::AngularCoordinate::Survey:
      lon_label = "Eta";
      lat_label = "Lambda";
      break;
    case Stomp::AngularCoordinate::Equatorial:
      lon_label = "RA";
      lat_label = "DEC";
      break;
    case Stomp::AngularCoordinate::Galactic:
      lon_label = "l";
      lat_label = "b";
      break;
    }

    if (width() < height()) {
      minor_tick_length_ = static_cast<int>(width()/100);
    } else {
      minor_tick_length_ = static_cast<int>(height()/100);
    }
    major_tick_length_ = 2*minor_tick_length_;
    minor_ticks_per_major_ = 4;

    uint16_t lon_tick_idx = 0;
    uint16_t n_lon_ticks =
      static_cast<int>((lonmax_ - lonmin_)/tick_divisions_[lon_tick_idx]);
    while ((n_lon_ticks > 10) && (lon_tick_idx < n_tick_divisions_ - 1)) {
      lon_tick_idx++;
      n_lon_ticks =
	static_cast<int>((lonmax_ - lonmin_)/tick_divisions_[lon_tick_idx]);
    }
    double lon_tick_step =
      tick_divisions_[lon_tick_idx]/(minor_ticks_per_major_ + 1.0);

    tick_precision_ = 4;
    if (Stomp::DoubleGE(tick_divisions_[lon_tick_idx], 0.01))
      tick_precision_ = 3;
    if (Stomp::DoubleGE(tick_divisions_[lon_tick_idx], 0.1))
      tick_precision_ = 2;
    if (Stomp::DoubleGE(tick_divisions_[lon_tick_idx], 1.0))
      tick_precision_ = 1;
    if (Stomp::DoubleGE(tick_divisions_[lon_tick_idx], 10.0))
      tick_precision_ = 0;

    double lon = -180.0;
    uint16_t major_tick_counter = minor_ticks_per_major_;
    while (Stomp::DoubleLE(lon, 360.0)) {
      if (Stomp::DoubleGE(lon, lonmin_) &&
	  Stomp::DoubleLE(lon, lonmax_)) {
	if (major_tick_counter == minor_ticks_per_major_) {
	  if (show_coordinates_) {
	    _drawLonTick(&painter, lon, true);
	    if (!Stomp::DoubleEQ(lon, lonmin_) &&
		!Stomp::DoubleEQ(lon, lonmax_))
	      _drawLonTickLabel(&painter, lon);
	  }
	  if (show_grid_) _drawLonGrid(&painter, lon);
	} else {
	  if (show_coordinates_) _drawLonTick(&painter, lon, false);
	}
      }
      if (major_tick_counter == minor_ticks_per_major_) {
	major_tick_counter = 0;
      } else {
	major_tick_counter++;
      }
      lon += lon_tick_step;
    }
    if (show_coordinates_) _drawLonLabel(&painter, lon_label);

    uint16_t lat_tick_idx = 0;
    uint16_t n_lat_ticks =
      static_cast<int>((latmax_ - latmin_)/tick_divisions_[lat_tick_idx]);
    while ((n_lat_ticks > 10) && (lat_tick_idx < n_tick_divisions_ - 1)) {
      lat_tick_idx++;
      n_lat_ticks =
	static_cast<int>((latmax_ - latmin_)/tick_divisions_[lat_tick_idx]);
    }
    double lat_tick_step =
      tick_divisions_[lat_tick_idx]/(minor_ticks_per_major_ + 1.0);

    tick_precision_ = 4;
    if (Stomp::DoubleGE(tick_divisions_[lat_tick_idx], 0.01))
      tick_precision_ = 3;
    if (Stomp::DoubleGE(tick_divisions_[lat_tick_idx], 0.1))
      tick_precision_ = 2;
    if (Stomp::DoubleGE(tick_divisions_[lat_tick_idx], 1.0))
      tick_precision_ = 1;
    if (Stomp::DoubleGE(tick_divisions_[lat_tick_idx], 10.0))
      tick_precision_ = 0;

    double lat = -90.0;
    major_tick_counter = minor_ticks_per_major_;
    while (Stomp::DoubleLE(lat, 90.0)) {
      if (Stomp::DoubleGE(lat, latmin_) &&
	  Stomp::DoubleLE(lat, latmax_)) {
	if (major_tick_counter == minor_ticks_per_major_) {
	  if (show_coordinates_) {
	    _drawLatTick(&painter, lat, true);
	    if (!Stomp::DoubleEQ(lat, latmin_) &&
		!Stomp::DoubleEQ(lat, latmax_))
	      _drawLatTickLabel(&painter, lat);
	  }
	  if (show_grid_) _drawLatGrid(&painter, lat);
	} else {
	  if (show_coordinates_) _drawLatTick(&painter, lat, false);
	}
      }
      if (major_tick_counter == minor_ticks_per_major_) {
	major_tick_counter = 0;
      } else {
	major_tick_counter++;
      }
      lat += lat_tick_step;
    }
    if (show_coordinates_) _drawLatLabel(&painter, lat_label);

    accessed_device = painter.end();
  }

  return accessed_device;
}

void RenderArea::_drawLonTick(QPainter *painter, double longitude,
			      bool major_tick) {
  uint16_t tick_length = minor_tick_length_;
  if (major_tick) tick_length = major_tick_length_;

  QPointF lower_point = QPointF(_lonToX(longitude),
				static_cast<qreal>(height() - 1));
  QPointF upper_point = QPointF(_lonToX(longitude),
				static_cast<qreal>(height() - 1 - tick_length));
  painter->drawLine(QLineF(lower_point, upper_point));

  lower_point = QPointF(_lonToX(longitude), 0.0);
  upper_point = QPointF(_lonToX(longitude), static_cast<qreal>(tick_length));
  painter->drawLine(QLineF(lower_point, upper_point));
}

void RenderArea::_drawLonTickLabel(QPainter* painter, double longitude) {
  // Make the label box as high and wide as the tick mark.
  qreal label_width = static_cast<qreal>(3*major_tick_length_);
  qreal label_height = static_cast<qreal>(major_tick_length_);

  // Center the box above the tick with a 1 pixel gap between the two.
  QPointF upper_left =
    QPointF(_lonToX(longitude) - 0.5*label_width,
	    static_cast<qreal>(height() - 2 - 2*label_height));
  QPointF lower_right =
    QPointF(_lonToX(longitude) + 0.5*label_width,
	    static_cast<qreal>(height() - 2 - label_height));
  painter->drawText(QRectF(upper_left, lower_right),
		    Qt::AlignCenter | Qt::TextDontClip,
		    QString::number(longitude, 'f', tick_precision_));
}

void RenderArea::_drawLatTick(QPainter *painter, double latitude,
			      bool major_tick) {
  uint16_t tick_length = minor_tick_length_;
  if (major_tick) tick_length = major_tick_length_;

  QPointF left_point = QPointF(0.0, _latToY(latitude));
  QPointF right_point = QPointF(static_cast<qreal>(tick_length),
				_latToY(latitude));
  painter->drawLine(QLineF(left_point, right_point));

  left_point = QPointF(static_cast<qreal>(width() - 1), _latToY(latitude));
  right_point = QPointF(static_cast<qreal>(width()- 1 - tick_length),
			_latToY(latitude));
  painter->drawLine(QLineF(left_point, right_point));
}

void RenderArea::_drawLatTickLabel(QPainter *painter, double latitude) {
  // Make the label box as high and wide as the tick mark and 4 times.
  qreal label_width = static_cast<qreal>(3*major_tick_length_);
  qreal label_height = static_cast<qreal>(major_tick_length_);

  // In order to draw the text rotated properly, we need to do a bit of
  // coordinate system tweaking.  After this translation and rotation, the
  // origin should be one pixel above the tick mark.
  painter->save();
  painter->translate(2*major_tick_length_ + 1.0, _latToY(latitude));
  painter->rotate(90.0);

  // Now we draw the text in native coordinates
  QPointF upper_left = QPointF(-0.5*label_width, label_height);
  QPointF lower_right = QPointF(0.5*label_width, 0.0);
  painter->drawText(QRectF(upper_left, lower_right),
		    Qt::AlignCenter | Qt::TextDontClip,
		    QString::number(latitude, 'f', tick_precision_));

  // Now reset to the original coordinate system.
  painter->restore();
}

void RenderArea::_drawLonLabel(QPainter* painter, QString& label) {
  // Make the label box as high as the tick mark and 4 times as wide.
  qreal label_width = static_cast<qreal>(5.0*major_tick_length_);
  qreal label_height = static_cast<qreal>(2.0*major_tick_length_);

  // Center the box above the tick with a 1 pixel gap between the two.
  double longitude = 0.5*(lonmin_ + lonmax_);
  QPointF upper_left =
    QPointF(_lonToX(longitude) - 0.5*label_width,
	    static_cast<qreal>(height() - 2 -
			       label_height - 2*major_tick_length_));
  QPointF lower_right =
    QPointF(_lonToX(longitude) + 0.5*label_width,
	    static_cast<qreal>(height() - 2 - 2*major_tick_length_));
  painter->drawText(QRectF(upper_left, lower_right), Qt::AlignCenter, label);
}

void RenderArea::_drawLatLabel(QPainter* painter, QString& label) {
  // Make the label box 1.5x as high as the tick mark and 4 times as wide.
  qreal label_width = static_cast<qreal>(5.0*major_tick_length_);
  qreal label_height = static_cast<qreal>(2.0*major_tick_length_);

  // In order to draw the text rotated properly, we need to do a bit of
  // coordinate system tweaking.  After this translation and rotation, the
  // origin should be one pixel above the tick mark.
  double latitude = 0.5*(latmin_ + latmax_);
  painter->save();
  painter->translate(4.0*major_tick_length_ + 2.0, _latToY(latitude));
  painter->rotate(90.0);

  // Now we draw the text in native coordinates
  QPointF upper_left = QPointF(-0.5*label_width, label_height);
  QPointF lower_right = QPointF(0.5*label_width, 0.0);
  painter->drawText(QRectF(upper_left, lower_right),
		    Qt::AlignCenter | Qt::TextDontClip, label);

  // Now reset to the original coordinate system.
  painter->restore();
}

void RenderArea::_drawLonGrid(QPainter* painter, double longitude) {
  QPolygonF poly_line;

  double input_longitude = longitude;

  uint16_t n_steps = 50;
  double lat_step = (latmax_ - latmin_)/n_steps;
  for (uint16_t i=0;i<=n_steps;i++){
    double latitude = latmin_ + i*lat_step;
    _enforceBounds(longitude, latitude);
    if (aitoff_projection_) _cartesianToAitoff(longitude, latitude);
    poly_line << QPointF(_lonToX(longitude), _latToY(latitude));
    longitude = input_longitude;
  }

  QPen pen = QPen(Qt::DashLine);
  pen.setColor(Qt::darkGray);
  painter->setPen(pen);
  painter->drawPolyline(poly_line);
  painter->setPen(QPen(Qt::darkGray));
}

void RenderArea::_drawLatGrid(QPainter* painter, double latitude) {
  QPolygonF poly_line;

  double input_latitude = latitude;

  uint16_t n_steps = 50;
  double lon_step = (lonmax_ - lonmin_)/n_steps;
  for (uint16_t i=0;i<=n_steps;i++){
    double longitude = lonmin_ + i*lon_step;
    _enforceBounds(longitude, latitude);
    if (aitoff_projection_) _cartesianToAitoff(longitude, latitude);
    poly_line << QPointF(_lonToX(longitude), _latToY(latitude));
    latitude = input_latitude;
  }

  QPen pen = QPen(Qt::DashLine);
  pen.setColor(Qt::darkGray);
  painter->setPen(pen);
  painter->drawPolyline(poly_line);
  painter->setPen(QPen(Qt::darkGray));
}

void RenderArea::_findNewWeightBounds(Stomp::Map* stomp_map) {
  setWeightRange(stomp_map->MinWeight(), stomp_map->MaxWeight());
}

void RenderArea::_findNewWeightBounds(Stomp::WAngularVector& ang) {
  double weight_min = ang[0].Weight();
  double weight_max = ang[0].Weight();

  for (uint32_t i=1;i<ang.size();i++) {
    if (weight_min > ang[i].Weight()) weight_min = ang[i].Weight();
    if (weight_max < ang[i].Weight()) weight_max = ang[i].Weight();
  }

  setPointsWeightRange(weight_min, weight_max);
}

void RenderArea::_findNewImageBounds(Stomp::Map* stomp_map) {
  if (!stomp_map->Empty()) {
    lonmax_ = -200.0;
    lonmin_ = 400.0;
    latmax_ = -200.0;
    latmin_ = 200.0;

    Stomp::PixelVector pix;
    if (stomp_map->Area() < Stomp::HPixArea) {
      // If the Map area is less than that of a single superpixel, then we
      // probably have a small enough map that we can check all of the pixels
      // directly and find the exact bounds.
      stomp_map->Pixels(pix);
    } else {
      // Otherwise, we'll get the Coverage pixels for the Map and use those
      // to work out our bounds.
      stomp_map->Coverage(pix);
    }

    // quick check against the possibility that we've got a full-sky map
    if (stomp_map->Area() > 0.2*4.0*Stomp::Pi*Stomp::StradToDeg) {
      setFullSky(true);
      setFullSky(false);
    } else {
      for (Stomp::PixelIterator iter=pix.begin();iter!=pix.end();++iter) {
	double lon = 0.0, lat = 0.0;
	double buffer = sqrt(iter->Area());
	switch (sphere_) {
	case Stomp::AngularCoordinate::Survey:
	  lon = iter->Eta();
	  lat = iter->Lambda();
	  break;
	case Stomp::AngularCoordinate::Equatorial:
	  lon = iter->RA();
	  lat = iter->DEC();
	  break;
	case Stomp::AngularCoordinate::Galactic:
	  lon = iter->GalLon();
	  lat = iter->GalLat();
	  break;
	}
	if (lon + buffer > lonmax_) lonmax_ = lon + buffer;
	if (lon - buffer < lonmin_) lonmin_ = lon - buffer;
	if (lat + buffer > latmax_) latmax_ = lat + buffer;
	if (lat - buffer < latmin_) latmin_ = lat - buffer;
      }
      _enforceBounds(lonmax_, latmax_);
      _enforceBounds(lonmin_, latmin_);
    }
  } else {
    setFullSky(true);
    setFullSky(false);
  }
}

void RenderArea::_findNewImageBounds(Stomp::WAngularVector& ang) {
  double latitude, longitude;

  switch (sphere_) {
  case Stomp::AngularCoordinate::Survey:
    longitude = ang[0].Eta();
    latitude = ang[0].Lambda();
    break;
  case Stomp::AngularCoordinate::Equatorial:
    longitude = ang[0].RA();
    latitude = ang[0].DEC();
    break;
  case Stomp::AngularCoordinate::Galactic:
    longitude = ang[0].GalLon();
    latitude = ang[0].GalLat();
    break;
  }

  lonmin_ = longitude;
  lonmax_ = longitude;
  latmin_ = latitude;
  latmax_ = latitude;

  for (uint32_t i=1;i<ang.size();i++) {
    switch (sphere_) {
    case Stomp::AngularCoordinate::Survey:
      longitude = ang[i].Eta();
      latitude = ang[i].Lambda();
      break;
    case Stomp::AngularCoordinate::Equatorial:
      longitude = ang[i].RA();
      latitude = ang[i].DEC();
      break;
    case Stomp::AngularCoordinate::Galactic:
      longitude = ang[i].GalLon();
      latitude = ang[i].GalLat();
      break;
    }

    if (lonmin_ > longitude) lonmin_ = longitude;
    if (lonmax_ < longitude) lonmax_ = longitude;
    if (latmin_ > latitude) latmin_ = latitude;
    if (latmax_ < latitude) latmax_ = latitude;
  }
}

void RenderArea::_findNewMaxResolution(Stomp::Map* stomp_map) {
  setMaxResolution(stomp_map->MaxResolution());
}

void RenderArea::_checkBounds() {
  continuous_longitude_ = true;
  if (Stomp::DoubleLT(latmin_, -90.0)) latmin_ = -90.0;
  if (Stomp::DoubleGT(latmax_, 90.0)) latmax_ = 90.0;
  if (sphere_ == Stomp::AngularCoordinate::Survey) {
    if (Stomp::DoubleLT(lonmin_, -180.0) ||
	Stomp::DoubleGT(lonmax_, 180.0)) continuous_longitude_ = false;
    while (Stomp::DoubleLT(lonmin_, -180.0)) lonmin_ += 360.0;
    while (Stomp::DoubleGT(lonmax_, 180.0)) lonmax_ -= 360.0;
  } else {
    if (Stomp::DoubleLT(lonmin_, 0.0) ||
	Stomp::DoubleGT(lonmax_, 360.0)) continuous_longitude_ = false;
    while (Stomp::DoubleLT(lonmin_, 0.0)) lonmin_ += 360.0;
    while (Stomp::DoubleGT(lonmax_, 360.0)) lonmax_ -= 360.0;
  }
  if (Stomp::DoubleGT(lonmin_, lonmax_)) continuous_longitude_ = false;
}

void RenderArea::_enforceBounds(double& lon, double& lat) {
  if (Stomp::DoubleGE(lat, 90.0)) lat = 90.0 - 0.0000001;
  if (Stomp::DoubleLE(lat, -90.0)) lat = -90.0 + 0.0000001;

  if (sphere_ == Stomp::AngularCoordinate::Survey) {
    if (Stomp::DoubleGE(lon, 180.0)) lon = 180.0 - 0.0000001;
    if (Stomp::DoubleLE(lon, -180.0)) lon = -180.0 + 0.0000001;
  } else {
    if (Stomp::DoubleGE(lon, 360.0)) lon = 360.0 - 0.0000001;
    if (Stomp::DoubleLE(lon, 0.0)) lon = 0.0 + 0.0000001;
  }
}

void RenderArea::_pixelToPoint(uint32_t pixel_x, uint32_t pixel_y,
			       uint32_t resolution, QPointF& point) {
  uint32_t nx = resolution*Stomp::Nx0;
  uint32_t ny = resolution*Stomp::Ny0;

  double lat_center =
    90.0 - Stomp::RadToDeg*acos(1.0-2.0*(pixel_y+0.5)/ny);
  double lon_center =
    Stomp::RadToDeg*(2.0*Stomp::Pi*(pixel_x+0.5))/nx +
    Stomp::EtaOffSet;
  if (Stomp::DoubleGT(lon_center, 180.0)) lon_center -= 360.0;
  if (Stomp::DoubleLT(lon_center, -180.0)) lon_center += 360.0;

  double tmp_lon, tmp_lat;

  switch (sphere_) {
  case Stomp::AngularCoordinate::Survey:
    break;
  case Stomp::AngularCoordinate::Equatorial:
    Stomp::AngularCoordinate::SurveyToEquatorial(lat_center, lon_center,
						 tmp_lon, tmp_lat);
    lon_center = tmp_lon;
    lat_center = tmp_lat;
    break;
  case Stomp::AngularCoordinate::Galactic:
    Stomp::AngularCoordinate::SurveyToGalactic(lat_center, lon_center,
					       tmp_lon, tmp_lat);
    lon_center = tmp_lon;
    lat_center = tmp_lat;
    break;
  }

  _enforceBounds(lon_center, lat_center);

  if (aitoff_projection_) _cartesianToAitoff(lon_center, lat_center);

  point.setX(_lonToX(lon_center));
  point.setY(_latToY(lat_center));
}

void RenderArea::_pixelToPolygon(uint32_t pixel_x, uint32_t pixel_y,
				 uint32_t resolution, QPolygonF& polygon) {
  if (!polygon.empty()) polygon.clear();

  uint32_t ratio = max_resolution_/resolution;

  uint32_t nx = resolution*Stomp::Nx0;
  uint32_t ny = resolution*Stomp::Ny0;

  double lat_center =
    90.0 - Stomp::RadToDeg*acos(1.0-2.0*(pixel_y+0.5)/ny);
  double lon_center =
    Stomp::RadToDeg*(2.0*Stomp::Pi*(pixel_x+0.5))/nx +
    Stomp::EtaOffSet;
  if (Stomp::DoubleGT(lon_center, 180.0)) lon_center -= 360.0;
  if (Stomp::DoubleLT(lon_center, -180.0)) lon_center += 360.0;

  double min_lon = 0.00000001;
  double mid_lon = 180.0;
  double max_lon = 359.99999999;
  double tmp_lon, tmp_lat;

  switch (sphere_) {
  case Stomp::AngularCoordinate::Survey:
    min_lon = -179.99999999;
    mid_lon = 0.0;
    max_lon = 179.99999999;
    break;
  case Stomp::AngularCoordinate::Equatorial:
    Stomp::AngularCoordinate::SurveyToEquatorial(lat_center, lon_center,
						 tmp_lon, tmp_lat);
    lon_center = tmp_lon;
    lat_center = tmp_lat;
    break;
  case Stomp::AngularCoordinate::Galactic:
    Stomp::AngularCoordinate::SurveyToGalactic(lat_center, lon_center,
					       tmp_lon, tmp_lat);
    lon_center = tmp_lon;
    lat_center = tmp_lat;
    break;
  }

  // pixel bottom edge
  for (uint32_t m=0;m<ratio;m++) {
    double lat = 90.0 - Stomp::RadToDeg*acos(1.0-2.0*(pixel_y+1.0)/ny);
    double lon =
      Stomp::RadToDeg*(2.0*Stomp::Pi*(pixel_x+1.0*m/ratio))/nx +
      Stomp::EtaOffSet;
    if (Stomp::DoubleGT(lon, 180.0000001)) lon -= 360.0;
    if (Stomp::DoubleLT(lon, -180.0000001)) lon += 360.0;

    if (sphere_ == Stomp::AngularCoordinate::Equatorial) {
      Stomp::AngularCoordinate::SurveyToEquatorial(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }
    if (sphere_ == Stomp::AngularCoordinate::Galactic) {
      Stomp::AngularCoordinate::SurveyToGalactic(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }

    if ((lon_center > max_lon - 2.0) && (lon < mid_lon)) lon = max_lon;
    if ((lon_center < min_lon + 2.0) && (lon > mid_lon)) lon = min_lon;
    _enforceBounds(lon, lat);

    if (aitoff_projection_) _cartesianToAitoff(lon, lat);
    polygon << QPointF(_lonToX(lon), _latToY(lat));
  }

  // pixel right edge
  for (uint32_t m=0;m<ratio;m++) {
    double lat = 90.0 -
      Stomp::RadToDeg*acos(1.0-2.0*(pixel_y+1.0*(ratio-m)/ratio)/ny);
    double lon =
      Stomp::RadToDeg*(2.0*Stomp::Pi*(pixel_x+1.0))/nx +
      Stomp::EtaOffSet;
    if (Stomp::DoubleGT(lon, 180.0000001)) lon -= 360.0;
    if (Stomp::DoubleLT(lon, -180.0000001)) lon += 360.0;

    if (sphere_ == Stomp::AngularCoordinate::Equatorial) {
      Stomp::AngularCoordinate::SurveyToEquatorial(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }
    if (sphere_ == Stomp::AngularCoordinate::Galactic) {
      Stomp::AngularCoordinate::SurveyToGalactic(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }

    if ((lon_center > max_lon - 2.0) && (lon < mid_lon)) lon = max_lon;
    if ((lon_center < min_lon + 2.0) && (lon > mid_lon)) lon = min_lon;
    _enforceBounds(lon, lat);

    if (aitoff_projection_) _cartesianToAitoff(lon, lat);
    polygon << QPointF(_lonToX(lon), _latToY(lat));
  }

  // pixel top edge
  for (uint32_t m=0;m<ratio;m++) {
    double lat =
      90.0 - Stomp::RadToDeg*acos(1.0-2.0*(pixel_y+0.0)/ny);
    double lon = Stomp::RadToDeg*
      (2.0*Stomp::Pi*(pixel_x+1.0*(ratio-m)/ratio))/nx +
      Stomp::EtaOffSet;
    if (Stomp::DoubleGT(lon, 180.0000001)) lon -= 360.0;
    if (Stomp::DoubleLT(lon, -180.0000001)) lon += 360.0;

    if (sphere_ == Stomp::AngularCoordinate::Equatorial) {
      Stomp::AngularCoordinate::SurveyToEquatorial(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }
    if (sphere_ == Stomp::AngularCoordinate::Galactic) {
      Stomp::AngularCoordinate::SurveyToGalactic(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }

    if ((lon_center > max_lon - 2.0) && (lon < mid_lon)) lon = max_lon;
    if ((lon_center < min_lon + 2.0) && (lon > mid_lon)) lon = min_lon;
    _enforceBounds(lon, lat);

    if (aitoff_projection_) _cartesianToAitoff(lon, lat);
    polygon << QPointF(_lonToX(lon), _latToY(lat));
  }

  // pixel left edge
  for (uint32_t m=0;m<ratio;m++) {
    double lat =
      90.0 - Stomp::RadToDeg*acos(1.0-2.0*(pixel_y+1.0*m/ratio)/ny);
    double lon =
      Stomp::RadToDeg*(2.0*Stomp::Pi*(pixel_x+0.0))/nx +
      Stomp::EtaOffSet;
    if (Stomp::DoubleGT(lon, 180.0000001)) lon -= 360.0;
    if (Stomp::DoubleLT(lon, -180.0000001)) lon += 360.0;

    if (sphere_ == Stomp::AngularCoordinate::Equatorial) {
      Stomp::AngularCoordinate::SurveyToEquatorial(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }
    if (sphere_ == Stomp::AngularCoordinate::Galactic) {
      Stomp::AngularCoordinate::SurveyToGalactic(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }

    if ((lon_center > max_lon - 2.0) && (lon < mid_lon)) lon = max_lon;
    if ((lon_center < min_lon + 2.0) && (lon > mid_lon)) lon = min_lon;
    _enforceBounds(lon, lat);

    if (aitoff_projection_) _cartesianToAitoff(lon, lat);
    polygon << QPointF(_lonToX(lon), _latToY(lat));
  }

  // repeat first point to close the polygon
  double lat =
    90.0 - Stomp::RadToDeg*acos(1.0-2.0*(pixel_y+1.0)/ny);
  double lon =
    Stomp::RadToDeg*(2.0*Stomp::Pi*(pixel_x+0.0))/nx +
    Stomp::EtaOffSet;
  if (Stomp::DoubleGT(lon, 180.0000001)) lon -= 360.0;
  if (Stomp::DoubleLT(lon, -180.0000001)) lon += 360.0;

  if (sphere_ == Stomp::AngularCoordinate::Equatorial) {
    Stomp::AngularCoordinate::SurveyToEquatorial(lat, lon, tmp_lon, tmp_lat);
    lon = tmp_lon;
    lat = tmp_lat;
  }
  if (sphere_ == Stomp::AngularCoordinate::Galactic) {
    Stomp::AngularCoordinate::SurveyToGalactic(lat, lon, tmp_lon, tmp_lat);
    lon = tmp_lon;
    lat = tmp_lat;
  }

  if ((lon_center > max_lon - 2.0) && (lon < mid_lon)) lon = max_lon;
  if ((lon_center < min_lon + 2.0) && (lon > mid_lon)) lon = min_lon;
  _enforceBounds(lon, lat);

  if (aitoff_projection_) _cartesianToAitoff(lon, lat);
  polygon << QPointF(_lonToX(lon), _latToY(lat));
}

void RenderArea::_splitPixelToPolygons(uint32_t pixel_x, uint32_t pixel_y,
				       uint32_t resolution,
				       QPolygonF& left_polygon,
				       QPolygonF& right_polygon) {
  if (!left_polygon.empty()) left_polygon.clear();
  if (!right_polygon.empty()) right_polygon.clear();

  uint32_t ratio = max_resolution_/resolution;

  uint32_t nx = resolution*Stomp::Nx0;
  uint32_t ny = resolution*Stomp::Ny0;

  double lat_center =
    90.0 - Stomp::RadToDeg*acos(1.0-2.0*(pixel_y+0.5)/ny);
  double lon_center =
    Stomp::RadToDeg*(2.0*Stomp::Pi*(pixel_x+0.5))/nx +
    Stomp::EtaOffSet;
  if (Stomp::DoubleGT(lon_center, 180.0)) lon_center -= 360.0;
  if (Stomp::DoubleLT(lon_center, -180.0)) lon_center += 360.0;

  double min_lon = 0.000001;
  double mid_lon = 180.0;
  double max_lon = 359.999999;
  double tmp_lon, tmp_lat;

  switch (sphere_) {
  case Stomp::AngularCoordinate::Survey:
    min_lon = -179.99999999;
    mid_lon = 0.0;
    max_lon = 179.99999999;
    break;
  case Stomp::AngularCoordinate::Equatorial:
    Stomp::AngularCoordinate::SurveyToEquatorial(lat_center, lon_center,
						 tmp_lon, tmp_lat);
    lon_center = tmp_lon;
    lat_center = tmp_lat;
    break;
  case Stomp::AngularCoordinate::Galactic:
    Stomp::AngularCoordinate::SurveyToGalactic(lat_center, lon_center,
					       tmp_lon, tmp_lat);
    lon_center = tmp_lon;
    lat_center = tmp_lat;
    break;
  }

  // pixel bottom edge
  for (uint32_t m=0;m<ratio;m++) {
    double lat = 90.0 - Stomp::RadToDeg*acos(1.0-2.0*(pixel_y+1.0)/ny);
    double lon =
      Stomp::RadToDeg*(2.0*Stomp::Pi*(pixel_x+1.0*m/ratio))/nx +
      Stomp::EtaOffSet;
    if (Stomp::DoubleGT(lon, 180.0)) lon -= 360.0;
    if (Stomp::DoubleLT(lon, -180.0)) lon += 360.0;

    if (sphere_ == Stomp::AngularCoordinate::Equatorial) {
      Stomp::AngularCoordinate::SurveyToEquatorial(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }
    if (sphere_ == Stomp::AngularCoordinate::Galactic) {
      Stomp::AngularCoordinate::SurveyToGalactic(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }

    double left_lon = lon;
    double right_lon = lon;
    if (lon > mid_lon) {
      left_lon = min_lon;
    } else {
      right_lon = max_lon;
    }
    _enforceBounds(left_lon, lat);
    _enforceBounds(right_lon, lat);
    if (aitoff_projection_) {
      _cartesianToAitoff(left_lon, lat);
      _cartesianToAitoff(right_lon, lat);
    }

    left_polygon << QPointF(_lonToX(left_lon), _latToY(lat));
    right_polygon << QPointF(_lonToX(right_lon), _latToY(lat));
  }

  // pixel right edge
  for (uint32_t m=0;m<ratio;m++) {
    double lat = 90.0 -
      Stomp::RadToDeg*acos(1.0-2.0*(pixel_y+1.0*(ratio-m)/ratio)/ny);
    double lon =
      Stomp::RadToDeg*(2.0*Stomp::Pi*(pixel_x+1.0))/nx +
      Stomp::EtaOffSet;
    if (Stomp::DoubleGT(lon, 180.0)) lon -= 360.0;
    if (Stomp::DoubleLT(lon, -180.0)) lon += 360.0;

    if (sphere_ == Stomp::AngularCoordinate::Equatorial) {
      Stomp::AngularCoordinate::SurveyToEquatorial(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }
    if (sphere_ == Stomp::AngularCoordinate::Galactic) {
      Stomp::AngularCoordinate::SurveyToGalactic(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }

    double left_lon = lon;
    double right_lon = lon;
    if (lon > mid_lon) {
      left_lon = min_lon;
    } else {
      right_lon = max_lon;
    }
    _enforceBounds(left_lon, lat);
    _enforceBounds(right_lon, lat);
    if (aitoff_projection_) {
      _cartesianToAitoff(left_lon, lat);
      _cartesianToAitoff(right_lon, lat);
    }

    left_polygon << QPointF(_lonToX(left_lon), _latToY(lat));
    right_polygon << QPointF(_lonToX(right_lon), _latToY(lat));
  }

  // pixel top edge
  for (uint32_t m=0;m<ratio;m++) {
    double lat =
      90.0 - Stomp::RadToDeg*acos(1.0-2.0*(pixel_y+0.0)/ny);
    double lon = Stomp::RadToDeg*
      (2.0*Stomp::Pi*(pixel_x+1.0*(ratio-m)/ratio))/nx +
      Stomp::EtaOffSet;
    if (Stomp::DoubleGT(lon, 180.0)) lon -= 360.0;
    if (Stomp::DoubleLT(lon, -180.0)) lon += 360.0;

    if (sphere_ == Stomp::AngularCoordinate::Equatorial) {
      Stomp::AngularCoordinate::SurveyToEquatorial(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }
    if (sphere_ == Stomp::AngularCoordinate::Galactic) {
      Stomp::AngularCoordinate::SurveyToGalactic(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }

    double left_lon = lon;
    double right_lon = lon;
    if (lon > mid_lon) {
      left_lon = min_lon;
    } else {
      right_lon = max_lon;
    }
    _enforceBounds(left_lon, lat);
    _enforceBounds(right_lon, lat);
    if (aitoff_projection_) {
      _cartesianToAitoff(left_lon, lat);
      _cartesianToAitoff(right_lon, lat);
    }

    left_polygon << QPointF(_lonToX(left_lon), _latToY(lat));
    right_polygon << QPointF(_lonToX(right_lon), _latToY(lat));
  }

  // pixel left edge
  for (uint32_t m=0;m<ratio;m++) {
    double lat =
      90.0 - Stomp::RadToDeg*acos(1.0-2.0*(pixel_y+1.0*m/ratio)/ny);
    double lon =
      Stomp::RadToDeg*(2.0*Stomp::Pi*(pixel_x+0.0))/nx +
      Stomp::EtaOffSet;
    if (Stomp::DoubleGT(lon, 180.0)) lon -= 360.0;
    if (Stomp::DoubleLT(lon, -180.0)) lon += 360.0;

    if (sphere_ == Stomp::AngularCoordinate::Equatorial) {
      Stomp::AngularCoordinate::SurveyToEquatorial(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }
    if (sphere_ == Stomp::AngularCoordinate::Galactic) {
      Stomp::AngularCoordinate::SurveyToGalactic(lat, lon, tmp_lon, tmp_lat);
      lon = tmp_lon;
      lat = tmp_lat;
    }

    double left_lon = lon;
    double right_lon = lon;
    if (lon > mid_lon) {
      left_lon = min_lon;
    } else {
      right_lon = max_lon;
    }
    _enforceBounds(left_lon, lat);
    _enforceBounds(right_lon, lat);
    if (aitoff_projection_) {
      _cartesianToAitoff(left_lon, lat);
      _cartesianToAitoff(right_lon, lat);
    }

    left_polygon << QPointF(_lonToX(left_lon), _latToY(lat));
    right_polygon << QPointF(_lonToX(right_lon), _latToY(lat));
  }

  // repeat first point to close the polygon
  double lat =
    90.0 - Stomp::RadToDeg*acos(1.0-2.0*(pixel_y+1.0)/ny);
  double lon =
    Stomp::RadToDeg*(2.0*Stomp::Pi*(pixel_x+0.0))/nx +
    Stomp::EtaOffSet;
  if (Stomp::DoubleGT(lon, 180.0)) lon -= 360.0;
  if (Stomp::DoubleLT(lon, -180.0)) lon += 360.0;

  if (sphere_ == Stomp::AngularCoordinate::Equatorial) {
    Stomp::AngularCoordinate::SurveyToEquatorial(lat, lon, tmp_lon, tmp_lat);
    lon = tmp_lon;
    lat = tmp_lat;
  }
  if (sphere_ == Stomp::AngularCoordinate::Galactic) {
    Stomp::AngularCoordinate::SurveyToGalactic(lat, lon, tmp_lon, tmp_lat);
    lon = tmp_lon;
    lat = tmp_lat;
  }

  double left_lon = lon;
  double right_lon = lon;
  if (lon > mid_lon) {
    left_lon = min_lon;
  } else {
    right_lon = max_lon;
  }
  _enforceBounds(left_lon, lat);
  _enforceBounds(right_lon, lat);
  if (aitoff_projection_) {
    _cartesianToAitoff(left_lon, lat);
    _cartesianToAitoff(right_lon, lat);
  }

  left_polygon << QPointF(_lonToX(left_lon), _latToY(lat));
  right_polygon << QPointF(_lonToX(right_lon), _latToY(lat));
}

bool RenderArea::_angToPoint(Stomp::WeightedAngularCoordinate& ang,
			     QPointF& point) {
  double latitude, longitude;

  switch (sphere_) {
  case Stomp::AngularCoordinate::Survey:
    longitude = ang.Eta();
    latitude = ang.Lambda();
    break;
  case Stomp::AngularCoordinate::Equatorial:
    longitude = ang.RA();
    latitude = ang.DEC();
    break;
  case Stomp::AngularCoordinate::Galactic:
    longitude = ang.GalLon();
    latitude = ang.GalLat();
    break;
  }

  bool inside_bounds = false;

  if (Stomp::DoubleLE(longitude, lonmax_) &&
      Stomp::DoubleGE(longitude, lonmin_) &&
      Stomp::DoubleLE(latitude, latmax_) &&
      Stomp::DoubleGE(latitude, latmin_)) {
    inside_bounds = true;

    if (aitoff_projection_) _cartesianToAitoff(longitude, latitude);

    point.setX(_lonToX(longitude));
    point.setY(_latToY(latitude));
  }

  return inside_bounds;
}

qreal RenderArea::_lonToX(double longitude) {
  double n_pixel = static_cast<double>(width() - buffer_left_ - buffer_right_);
  double range = longitude - lonmin_;

  if (!continuous_longitude_) {
    range = 360.0 - lonmin_ + longitude;
    if (sphere_ == Stomp::AngularCoordinate::Survey) {
      if (Stomp::DoubleLT(longitude, lonmin_))
	range = lonmin_ - longitude;
    } else {
      if (Stomp::DoubleGT(longitude, lonmin_))
	range = longitude - lonmin_;
    }
  }

  double pixel_x = range*n_pixel/_lonRange();
  if (Stomp::DoubleLT(pixel_x, 0.0)) pixel_x = 0.0;
  if (Stomp::DoubleGE(pixel_x, n_pixel))
    pixel_x = n_pixel;

  return static_cast<qreal>(pixel_x + buffer_left_);
}

double RenderArea::_xToLon(qreal pixel_x) {
  double n_pixel = static_cast<double>(width() - buffer_left_ - buffer_right_);
  pixel_x -= static_cast<qreal>(buffer_left_);

  double range = static_cast<double>(_lonRange()*(pixel_x/n_pixel));

  double longitude = range + lonmin_;
  if (!continuous_longitude_) {
    longitude = lonmax_ + range;
    if (sphere_ == Stomp::AngularCoordinate::Survey) {
      if (Stomp::DoubleGE(longitude, 180.0)) longitude -= 360.0;
    } else {
      if (Stomp::DoubleGE(longitude, 360.0)) longitude -= 360.0;
    }
  }

  return longitude;
}

qreal RenderArea::_latToY(double latitude) {
  double n_pixel = static_cast<double>(height() - buffer_top_ - buffer_bottom_);
  double range = latmax_ - latitude;

  double pixel_y = range*n_pixel/_latRange();
  if (Stomp::DoubleLT(pixel_y, 0.0)) pixel_y = 0.0;
  if (Stomp::DoubleGE(pixel_y, n_pixel)) pixel_y = n_pixel;

  return static_cast<qreal>(pixel_y + buffer_top_);
}

double RenderArea::_yToLat(qreal pixel_y) {
  double n_pixel = static_cast<double>(height() - buffer_top_ - buffer_bottom_);
  pixel_y -= static_cast<qreal>(buffer_top_);

  double range = static_cast<double>(_latRange()*(pixel_y/n_pixel));
  double latitude = latmax_ - range;

  if (Stomp::DoubleLE(latitude, latmin_)) latitude = latmin_;
  if (Stomp::DoubleGE(latitude, latmax_)) latitude = latmax_;

  return latitude;
}

double RenderArea::_lonRange() {
  return (continuous_longitude_ ? lonmax_ - lonmin_ :
	  360.0 - lonmin_ + lonmax_);
}

double RenderArea::_latRange() {
  return latmax_ - latmin_;
}

double RenderArea::_lonCenter() {
  double lon_center = 0.5*(lonmin_ + lonmax_);
  if (!continuous_longitude_) {
    lon_center = lonmax_ + 0.5*_lonRange();
    if (sphere_ == Stomp::AngularCoordinate::Survey) {
      if (Stomp::DoubleGE(lon_center, 180.0)) lon_center -= 360.0;
    } else {
      if (Stomp::DoubleGE(lon_center, 360.0)) lon_center -= 360.0;
    }
  }

  return lon_center;
}

double RenderArea::_latCenter() {
  return 0.5*(latmin_ + latmax_);
}

void RenderArea::_cartesianToAitoff(double& longitude, double& latitude) {
  double sa = longitude;
  if ((sphere_ == Stomp::AngularCoordinate::Equatorial) ||
      (sphere_ == Stomp::AngularCoordinate::Galactic)) {
    /*    if (Stomp::DoubleLT(longitude, 180.0)) {
	  sa = longitude + 180.0;
	  } else {
	  sa = longitude - 180.0;
	  } */
    // if (Stomp::DoubleGE(sa, 180.0)) sa -= 360.0;
    sa -= 180.0;
  }
  double alpha2 = sa*Stomp::DegToRad/2.0;
  double delta = latitude*Stomp::DegToRad;
  double r2 = sqrt(2.0);
  double f = 2.0*r2/Stomp::Pi;

  double cdec = cos(delta);
  double denom = sqrt(1.0 + cdec*cos(alpha2))*(f*Stomp::DegToRad);

  longitude = cdec*sin(alpha2)*2.0*r2/denom;
  if ((sphere_ == Stomp::AngularCoordinate::Equatorial) ||
      (sphere_ == Stomp::AngularCoordinate::Galactic)) longitude += 180.0;
  latitude = sin(delta)*r2/denom;
}

double RenderArea::_normalizeWeight(double weight) {
  if (Stomp::DoubleGE(weight, weight_max_)) weight = weight_max_;
  if (Stomp::DoubleLE(weight, weight_min_)) weight = weight_min_;
  return (weight - weight_min_)/(weight_max_ - weight_min_);
}

double RenderArea::_normalizePointsWeight(double weight) {
  if (Stomp::DoubleGE(weight, points_weight_max_)) weight = points_weight_max_;
  if (Stomp::DoubleLE(weight, points_weight_min_)) weight = points_weight_min_;
  return
    (weight - points_weight_min_)/(points_weight_max_ - points_weight_min_);
}

double RenderArea::_renderPixelArea() {
  double n_pixel =
    static_cast<double>((width() - buffer_left_ - buffer_right_)*
			(height() - buffer_top_ - buffer_bottom_));

  double cartesian_area = _lonRange()*_latRange();

  return cartesian_area/n_pixel;
}

void RenderArea::_findRenderLevel(uint32_t target_pixels) {
  double cartesian_area = _lonRange()*_latRange();

  render_level_ = Stomp::HPixLevel;
  for (uint8_t level=Stomp::HPixLevel+1;level<=Stomp::MaxPixelLevel;level++) {
    uint32_t n_pixel = stomp_map_[level]->Size();

    if (cartesian_area < 10000.0) {
      n_pixel = 0;

      Stomp::Map* level_map = stomp_map_[level];
      for (Stomp::MapIterator iter=level_map->Begin();
	   iter!=level_map->End();level_map->Iterate(&iter)) {
	Stomp::PixelIterator pix = iter.second;      
	double latitude, longitude;

	switch (sphere_) {
	case Stomp::AngularCoordinate::Survey:
	  longitude = pix->Eta();
	  latitude = pix->Lambda();
	  break;
	case Stomp::AngularCoordinate::Equatorial:
	  longitude = pix->RA();
	  latitude = pix->DEC();
	  break;
	case Stomp::AngularCoordinate::Galactic:
	  longitude = pix->GalLon();
	  latitude = pix->GalLat();
	  break;
	}

	if (Stomp::DoubleLE(longitude, lonmax_) &&
	    Stomp::DoubleGE(longitude, lonmin_) &&
	    Stomp::DoubleLE(latitude, latmax_) &&
	    Stomp::DoubleGE(latitude, latmin_)) n_pixel++;
      }
    }
      
    if (n_pixel < target_pixels) {
      render_level_ = level;
    } else {
      break;
    }
  }
}

void RenderArea::_clearMaps() {
  delete stomp_map_[Stomp::MaxPixelLevel];
  for (uint8_t level=Stomp::MaxPixelLevel-1;level>=Stomp::HPixLevel;level--) {
    if (stomp_map_[level] != stomp_map_[level+1])
      delete stomp_map_[level];
  }

  stomp_map_.clear();
  good_map_ = false;
}

