#include <QtGui>
#include "render_area.h"

RenderArea::RenderArea(QWidget *parent) : QWidget(parent) {
  antialiased_ = false;
  auto_update_ = false;

  setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
  setMouseTracking(true);

  map_palette_.setWeightRange(0.0, 1.0);

  points_palette_.setWeightRange(0.0, 1.0);

  /*  Stomp::PixelVector pix;
      double dweight = 1.0/(4*Stomp::MaxSuperpixnum);
      for (uint32_t k=0;k<4*Stomp::MaxSuperpixnum;k++)
      pix.push_back(Stomp::Pixel(8, k, k*dweight));
      Stomp::Map* base_map = new Stomp::Map(pix);

      for (uint8_t level=Stomp::HPixLevel;level<=Stomp::MaxPixelLevel;level++)
      stomp_map_[level] = base_map;
  */
  good_map_ = false;
  good_points_ = false;

  fill_ = false;

  geom_.setImageBounds(0.0, 360.0, -90.0, 90.0,
		       Stomp::AngularCoordinate::Equatorial);
  geom_.setHeight(height());
  geom_.setWidth(width());

  setFullSky(false);
  useAitoffProjection(false);

  max_resolution_ = Stomp::HPixResolution;

  setBackgroundRole(QPalette::Base);
  setAutoFillBackground(true);

  background_color_ = Qt::white;
  map_pixmap_ = QPixmap(width(), height());
  map_pixmap_.fill(background_color_);

  qRegisterMetaType<QImage>("QImage");
  connect(&render_map_thread_, SIGNAL(renderedImage(QImage)),
	  this, SLOT(newMapImage(QImage)));
  connect(&render_map_thread_, SIGNAL(renderProgress(int)),
          this, SLOT(renderProgress(int)));
  connect(&render_points_thread_, SIGNAL(renderedImage(QImage)),
	  this, SLOT(newPointsImage(QImage)));
  connect(&render_points_thread_, SIGNAL(renderProgress(int)),
          this, SLOT(renderProgress(int)));

  connect(&read_map_thread_, SIGNAL(newBaseMap(Stomp::Map*)),
	  this, SLOT(newBaseMap(Stomp::Map*)));

  show_grid_ = false;
  show_coordinates_ = false;
  grid_pixmap_ = new QPixmap(width(), height());
  grid_pixmap_->fill(Qt::transparent);

  show_points_ = false;
  filter_points_ = false;
  points_pixmap_ = QPixmap(width(), height());
  points_pixmap_.fill(Qt::transparent);

  n_tick_divisions_ = 27;
  tick_divisions_ = new double[n_tick_divisions_];
  for (int i=0, scale=1;i<8;i++, scale*=10) {
    tick_divisions_[3*i + 0] = 0.0000001*scale;
    tick_divisions_[3*i + 1] = 0.0000002*scale;
    tick_divisions_[3*i + 2] = 0.0000005*scale;
  }
  tick_divisions_[24] = 10.0;
  tick_divisions_[25] = 15.0;
  tick_divisions_[26] = 30.0;

  major_tick_length_ = 0;
  minor_tick_length_ = 0;
  tick_precision_ = 4;
  minor_ticks_per_major_ = 4;

#ifndef QT_NO_CURSOR
  setCursor(Qt::CrossCursor);
#endif
}

QSize RenderArea::minimumSizeHint() const {
  return QSize(200, 100);
}

QSize RenderArea::sizeHint() const {
  return QSize(1024, 512);
}

double RenderArea::mapWeightMin() {
  return map_palette_.weightMin();
}

double RenderArea::mapWeightMax() {
  return map_palette_.weightMax();
}

double RenderArea::pointsWeightMin() {
  return points_palette_.weightMin();
}

double RenderArea::pointsWeightMax() {
  return points_palette_.weightMax();
}

double RenderArea::longitudeMin() {
  return geom_.longitudeMin();
}

double RenderArea::longitudeMax() {
  return geom_.longitudeMax();
}

double RenderArea::latitudeMin() {
  return geom_.latitudeMin();
}

double RenderArea::latitudeMax() {
  return geom_.latitudeMax();
}

uint32_t RenderArea::maxResolution() {
  return max_resolution_;
}

bool RenderArea::autoUpdating() {
  return auto_update_;
}

bool RenderArea::fullSky() {
  return geom_.fullSky();
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

double RenderArea::mouseLongitude() {
  return mouse_lon_;
}

double RenderArea::mouseLatitude() {
  return mouse_lat_;
}

void RenderArea::updatePixmap(bool send_update_call) {
  if (good_map_) {
    render_map_thread_.renderMap(base_map_, geom_, map_palette_,
				 max_resolution_, antialiased_, fill_);
  } else {
    map_pixmap_ = QPixmap(width(), height());
    map_pixmap_.fill(background_color_);
  }

  if (send_update_call)  updatePoints(false);
  updateGrid(false);
  if (send_update_call) update();
}

void RenderArea::newMapImage(const QImage& map_image) {
  map_pixmap_ = QPixmap::fromImage(map_image);

  update();
}

void RenderArea::updatePoints(bool send_update_call) {
  if (good_points_) {
    render_points_thread_.renderPoints(&ang_, geom_, points_palette_,
				       antialiased_);
  } else {
    points_pixmap_ = QPixmap(width(), height());
    points_pixmap_.fill(Qt::transparent);
  }

  if (send_update_call) update();
}

void RenderArea::newPointsImage(const QImage& points_image) {
  points_pixmap_ = QPixmap::fromImage(points_image);

  update();
}

void RenderArea::updateGrid(bool send_update_call) {
  delete grid_pixmap_;

  grid_pixmap_ = new QPixmap(size());
  grid_pixmap_->fill(Qt::transparent);

  paintGrid(grid_pixmap_);
  if (send_update_call) update();
}

bool RenderArea::readNewMap(QString& input_file) {
  std::string stl_input_file = std::string(input_file.toLatin1());
  bool file_exists = false;
  std::ifstream test_file(stl_input_file.c_str());

  if (good_map_) {
    base_map_->Clear();
    delete base_map_;
  }

  good_map_ = false;

  if (test_file) {
    file_exists = true;
    test_file.close();
    read_map_thread_.readMap(stl_input_file);
  }

  return file_exists;
}

void RenderArea::newBaseMap(Stomp::Map* base_map) {
  base_map_ = base_map;

  good_map_ = !base_map_->Empty();
  if (good_map_) {
    emit newStatus("Read base map.  Rendering...", 0);
    _findNewWeightBounds(base_map);
    _findNewImageBounds(base_map);
    _findNewMaxResolution(base_map);
    emit newMapParameters();
    updatePixmap();
  }
}

void RenderArea::renderProgress(int progress) {
  emit progressUpdate(progress);
}

bool RenderArea::readNewPointsFile(QString& input_file,
				   Stomp::AngularCoordinate::Sphere sphere,
				   bool weighted_points) {
  ang_.clear();

  std::string stl_input_file = std::string(input_file.toLatin1());

  int8_t weight_column = -1;
  if (weighted_points) weight_column = 2;

  good_points_ =
    Stomp::WeightedAngularCoordinate::ToWAngularVector(stl_input_file, ang_,
						       sphere, false, 0, 1,
						       weight_column);

  if (good_points_) {
    _findNewWeightBounds(ang_);
    if (!good_map_) _findNewImageBounds(ang_);
    emit newMapParameters();

    updatePoints();
  }

  return good_points_;
}

void RenderArea::clearPoints() {
  ang_.clear();
  updatePoints();
}

bool RenderArea::writeToPng(QString& output_file) {
  if (size() != map_pixmap_.size()) updatePixmap();

  QPixmap* output_pixmap = new QPixmap(width(), height());

  QPainter painter(output_pixmap);
  painter.drawPixmap(0, 0, map_pixmap_);
  if (show_points_) painter.drawPixmap(0, 0, points_pixmap_);
  if (show_grid_ || show_coordinates_) painter.drawPixmap(0, 0, *grid_pixmap_);

  bool painted_map = output_pixmap->save(output_file);

  return painted_map;
}

void RenderArea::setAntialiased(bool antialiased) {
  this->antialiased_ = antialiased;
  if (auto_update_) updatePixmap();
}

void RenderArea::setAutoUpdate(bool auto_update) {
  this->auto_update_ = auto_update;
}

void RenderArea::setMapWeightRange(double weight_min, double weight_max) {
  this->map_palette_.setWeightRange(weight_min, weight_max);
  if (auto_update_) updatePixmap();
}

void RenderArea::setPointsWeightRange(double weight_min, double weight_max) {
  this->points_palette_.setWeightRange(weight_min, weight_max);
  if (auto_update_) updatePixmap();
}

void RenderArea::setMapPalette(Palette::PaletteType palette_type) {
  this->map_palette_.initialize(palette_type);
  if (auto_update_) updatePixmap();
}

void RenderArea::setPointsPalette(Palette::PaletteType palette_type) {
  this->points_palette_.initialize(palette_type);
  if (auto_update_) updatePixmap();
}

void RenderArea::setVerticalBuffer(uint16_t buffer_pixels) {
  geom_.setVerticalBuffer(buffer_pixels);
  if (auto_update_) updatePixmap();
}

void RenderArea::setHorizontalBuffer(uint16_t buffer_pixels) {
  geom_.setHorizontalBuffer(buffer_pixels);
  if (auto_update_) updatePixmap();
}

void RenderArea::setTopBuffer(uint16_t buffer_pixels) {
  geom_.setTopBuffer(buffer_pixels);
  if (auto_update_) updatePixmap();
}

void RenderArea::setBottomBuffer(uint16_t buffer_pixels) {
  geom_.setBottomBuffer(buffer_pixels);
  if (auto_update_) updatePixmap();
}

void RenderArea::setLeftBuffer(uint16_t buffer_pixels) {
  geom_.setLeftBuffer(buffer_pixels);
  if (auto_update_) updatePixmap();
}

void RenderArea::setRightBuffer(uint16_t buffer_pixels) {
  geom_.setRightBuffer(buffer_pixels);
  if (auto_update_) updatePixmap();
}

void RenderArea::setImageBounds(double lonmin, double lonmax,
				double latmin, double latmax) {
  geom_.setImageBounds(lonmin, lonmax, latmin, latmax);
  if (auto_update_) updatePixmap();
}

void RenderArea::setImageBounds(double lonmin, double lonmax,
				double latmin, double latmax,
				Stomp::AngularCoordinate::Sphere sphere) {
  geom_.setImageBounds(lonmin, lonmax, latmin, latmax, sphere);
  if (auto_update_) updatePixmap();
}

void RenderArea::setFullSky(bool full_sky) {
  geom_.setFullSky(full_sky);
  if (full_sky && auto_update_) updatePixmap();
}

void RenderArea::setCoordinateSystem(Stomp::AngularCoordinate::Sphere sphere) {
  geom_.setCoordinateSystem(sphere);
  if (geom_.fullSky()) {
    setFullSky(true);
  } else {
    if (good_map_) {
      _findNewImageBounds(base_map_);
    } else {
      setFullSky(true);
      setFullSky(false);
    }
  }
  if (auto_update_) updatePixmap();
}

void RenderArea::zoomIn() {
  zoomIn(geom_.lonCenter(), geom_.latCenter());
}

void RenderArea::zoomOut() {
  zoomOut(geom_.lonCenter(), geom_.latCenter());
}

void RenderArea::zoomIn(double zoom_factor) {
  zoomIn(geom_.lonCenter(), geom_.latCenter(), zoom_factor);
}

void RenderArea::zoomOut(double zoom_factor) {
  zoomOut(geom_.lonCenter(), geom_.latCenter(), zoom_factor);
}

void RenderArea::fillPolygons(bool fill_polygons) {
  fill_ = fill_polygons;
  if (auto_update_) updatePixmap();
}

void RenderArea::useAitoffProjection(bool aitoff_projection) {
  geom_.useAitoffProjection(aitoff_projection);
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
  if (size() != map_pixmap_.size()) updatePixmap(false);
  if (size() != points_pixmap_.size()) updatePoints(false);

  QPainter painter(this);
  painter.drawPixmap(0, 0, map_pixmap_);
  if (show_points_) painter.drawPixmap(0, 0, points_pixmap_);
  if (show_grid_ || show_coordinates_) painter.drawPixmap(0, 0, *grid_pixmap_);
}

void RenderArea::mouseMoveEvent(QMouseEvent *event) {
  QPointF mouse_position = event->posF();

  if ((mouse_position.x() < width()) && (mouse_position.x() > 0) &&
      (mouse_position.y() < height()) && (mouse_position.y() > 0)) {
    mouse_lon_ = geom_.xToLon(mouse_position.x());
    mouse_lat_ = geom_.yToLat(mouse_position.y());
    emit newMousePosition();
  }
}

void RenderArea::mouseDoubleClickEvent(QMouseEvent *event) {
  QPointF new_center = event->posF();

  if (!geom_.aitoffProjection())
    zoomIn(geom_.xToLon(new_center.x()), geom_.yToLat(new_center.y()));
}

void RenderArea::resizeEvent(QResizeEvent *event) {
  geom_.setSize(event->size());
}

void RenderArea::zoomIn(double new_longitude_center,
			double new_latitude_center,
			double zoom_factor) {
  // After zooming in, the longitude range (2*lon_span), should decrease by
  // 1/zoom_factor.
  double lon_span = 0.5*geom_.lonRange()/zoom_factor;
  double lat_span = 0.5*geom_.latRange()/zoom_factor;

  geom_.setImageBounds(new_longitude_center - lon_span,
		       new_longitude_center + lon_span,
		       new_latitude_center - lat_span,
		       new_latitude_center + lat_span);

  updatePixmap();

  emit newMapParameters();
}

void RenderArea::zoomOut(double new_longitude_center,
			 double new_latitude_center,
			 double zoom_factor) {
  // After zooming out, the longitude range (2*lon_span), should increase by
  // zoom_factor.
  double lon_span = 0.5*geom_.lonRange()*zoom_factor;
  double lat_span = 0.5*geom_.latRange()*zoom_factor;

  if (Stomp::DoubleLE(lon_span, 180.0) &&
      Stomp::DoubleLE(lat_span, 90.0)) {
    // Provided that zooming out won't max out either of our coordinates, we
    // can safely increase each range.
    geom_.setImageBounds(new_longitude_center - lon_span,
			 new_longitude_center + lon_span,
			 new_latitude_center - lat_span,
			 new_latitude_center + lat_span);
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
      double lonmin = geom_.longitudeMin(), lonmax = geom_.longitudeMax();
      double latmin = geom_.latitudeMin(), latmax = geom_.latitudeMax();
      if (Stomp::DoubleLE(lon_span, 180.0)) {
	lonmin = new_longitude_center - lon_span;
	lonmax = new_longitude_center + lon_span;
      } else {
	if (geom_.sphere() == Stomp::AngularCoordinate::Survey) {
	  lonmin = -180.0;
	  lonmax = 180.0;
	} else {
	  lonmin = 0.0;
	  lonmax = 360.0;
	}
      }
      if (Stomp::DoubleLE(lat_span, 90.0)) {
	latmin = new_latitude_center - lat_span;
	latmax = new_latitude_center + lat_span;
      } else {
	latmin = -90.0;
	latmax = 90.0;
      }
      geom_.setImageBounds(lonmin, lonmax, latmin, latmax);
    }
  }

  updatePixmap();

  emit newMapParameters();
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
	if (geom_.angToPoint(ang_[i], point)) {
	  QBrush pixel_brush = QBrush(points_palette_.color(ang_[i].Weight()));
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

    switch (geom_.sphere()) {
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

    uint32_t lon_tick_idx = 0;
    uint32_t n_lon_ticks =
      static_cast<uint32_t>(geom_.lonRange()/tick_divisions_[lon_tick_idx]);
    while ((n_lon_ticks > 10) && (lon_tick_idx < n_tick_divisions_ - 1)) {
      lon_tick_idx++;
      n_lon_ticks =
	static_cast<uint32_t>(geom_.lonRange()/tick_divisions_[lon_tick_idx]);
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
    uint32_t major_tick_counter = minor_ticks_per_major_;
    while (Stomp::DoubleLE(lon, 360.0)) {
      if (Stomp::DoubleGE(lon, geom_.longitudeMin()) &&
	  Stomp::DoubleLE(lon, geom_.longitudeMax())) {
	if (major_tick_counter == minor_ticks_per_major_) {
	  if (show_coordinates_) {
	    _drawLonTick(&painter, lon, true);
	    if (!Stomp::DoubleEQ(lon, geom_.longitudeMin()) &&
		!Stomp::DoubleEQ(lon, geom_.longitudeMax()))
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

    uint32_t lat_tick_idx = 0;
    uint32_t n_lat_ticks =
      static_cast<uint32_t>(geom_.latRange()/tick_divisions_[lat_tick_idx]);
    while ((n_lat_ticks > 10) && (lat_tick_idx < n_tick_divisions_ - 1)) {
      lat_tick_idx++;
      n_lat_ticks =
	static_cast<uint32_t>(geom_.latRange()/tick_divisions_[lat_tick_idx]);
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
      if (Stomp::DoubleGE(lat, geom_.latitudeMin()) &&
	  Stomp::DoubleLE(lat, geom_.latitudeMax())) {
	if (major_tick_counter == minor_ticks_per_major_) {
	  if (show_coordinates_) {
	    _drawLatTick(&painter, lat, true);
	    if (!Stomp::DoubleEQ(lat, geom_.latitudeMin()) &&
		!Stomp::DoubleEQ(lat, geom_.latitudeMax()))
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

  QPointF lower_point = QPointF(geom_.lonToX(longitude),
				static_cast<qreal>(height() - 1));
  QPointF upper_point = QPointF(geom_.lonToX(longitude),
				static_cast<qreal>(height() - 1 - tick_length));
  painter->drawLine(QLineF(lower_point, upper_point));

  lower_point = QPointF(geom_.lonToX(longitude), 0.0);
  upper_point = QPointF(geom_.lonToX(longitude),
			static_cast<qreal>(tick_length));
  painter->drawLine(QLineF(lower_point, upper_point));
}

void RenderArea::_drawLonTickLabel(QPainter* painter, double longitude) {
  // Make the label box as high and wide as the tick mark.
  qreal label_width = static_cast<qreal>(3*major_tick_length_);
  qreal label_height = static_cast<qreal>(major_tick_length_);

  // Center the box above the tick with a 1 pixel gap between the two.
  QPointF upper_left =
    QPointF(geom_.lonToX(longitude) - 0.5*label_width,
	    static_cast<qreal>(height() - 2 - 2*label_height));
  QPointF lower_right =
    QPointF(geom_.lonToX(longitude) + 0.5*label_width,
	    static_cast<qreal>(height() - 2 - label_height));
  painter->drawText(QRectF(upper_left, lower_right),
		    Qt::AlignCenter | Qt::TextDontClip,
		    QString::number(longitude, 'f', tick_precision_));
}

void RenderArea::_drawLatTick(QPainter *painter, double latitude,
			      bool major_tick) {
  uint16_t tick_length = minor_tick_length_;
  if (major_tick) tick_length = major_tick_length_;

  QPointF left_point = QPointF(0.0, geom_.latToY(latitude));
  QPointF right_point = QPointF(static_cast<qreal>(tick_length),
				geom_.latToY(latitude));
  painter->drawLine(QLineF(left_point, right_point));

  left_point = QPointF(static_cast<qreal>(width() - 1), geom_.latToY(latitude));
  right_point = QPointF(static_cast<qreal>(width()- 1 - tick_length),
			geom_.latToY(latitude));
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
  painter->translate(2*major_tick_length_ + 1.0, geom_.latToY(latitude));
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
  double longitude = geom_.lonCenter();
  QPointF upper_left =
    QPointF(geom_.lonToX(longitude) - 0.5*label_width,
	    static_cast<qreal>(height() - 2 -
			       label_height - 2*major_tick_length_));
  QPointF lower_right =
    QPointF(geom_.lonToX(longitude) + 0.5*label_width,
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
  double latitude = geom_.latCenter();
  painter->save();
  painter->translate(4.0*major_tick_length_ + 2.0, geom_.latToY(latitude));
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
  double lat_step = geom_.latRange()/n_steps;
  for (uint16_t i=0;i<=n_steps;i++){
    double latitude = geom_.latitudeMin() + i*lat_step;
    geom_.enforceBounds(longitude, latitude);
    if (geom_.aitoffProjection()) geom_.cartesianToAitoff(longitude, latitude);
    poly_line << QPointF(geom_.lonToX(longitude), geom_.latToY(latitude));
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
  double lon_step = geom_.lonRange()/n_steps;
  for (uint16_t i=0;i<=n_steps;i++){
    double longitude = geom_.longitudeMin() + i*lon_step;
    geom_.enforceBounds(longitude, latitude);
    if (geom_.aitoffProjection()) geom_.cartesianToAitoff(longitude, latitude);
    poly_line << QPointF(geom_.lonToX(longitude), geom_.latToY(latitude));
    latitude = input_latitude;
  }

  QPen pen = QPen(Qt::DashLine);
  pen.setColor(Qt::darkGray);
  painter->setPen(pen);
  painter->drawPolyline(poly_line);
  painter->setPen(QPen(Qt::darkGray));
}

void RenderArea::_findNewWeightBounds(Stomp::Map* stomp_map) {
  setMapWeightRange(stomp_map->MinWeight(), stomp_map->MaxWeight());
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
    double lonmax = -200.0, lonmin = 400.0;
    double latmax = -200.0, latmin = 200.0;

    Stomp::PixelVector pix;
    if (stomp_map->Area() < Stomp::HPixArea) {
      // If the Map area is less than that of a single superpixel, then we
      // probably have a small enough map that we can check all of the pixels
      // directly and find the exact bounds.
      stomp_map->Pixels(pix);
    } else {
      // Otherwise, we'll get the Coverage pixels for the Map and use those
      // to work out our bounds.
      stomp_map->Coverage(pix, Stomp::HPixResolution, false);
    }

    // quick check against the possibility that we've got a full-sky map
    if (stomp_map->Area() > 0.2*4.0*Stomp::Pi*Stomp::StradToDeg) {
      setFullSky(true);
      setFullSky(false);
    } else {
      for (Stomp::PixelIterator iter=pix.begin();iter!=pix.end();++iter) {
	double lon = 0.0, lat = 0.0;
	double buffer = sqrt(iter->Area());
	switch (geom_.sphere()) {
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
	if (lon + buffer > lonmax) lonmax = lon + buffer;
	if (lon - buffer < lonmin) lonmin = lon - buffer;
	if (lat + buffer > latmax) latmax = lat + buffer;
	if (lat - buffer < latmin) latmin = lat - buffer;
      }
      geom_.enforceBounds(lonmax, latmax);
      geom_.enforceBounds(lonmin, latmin);
      geom_.setImageBounds(lonmin, lonmax, latmin, latmax);
    }
  } else {
    setFullSky(true);
    setFullSky(false);
  }
}

void RenderArea::_findNewImageBounds(Stomp::WAngularVector& ang) {
  double latitude, longitude;

  switch (geom_.sphere()) {
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

  double lonmin = longitude, lonmax = longitude;
  double latmin = latitude, latmax = latitude;

  for (uint32_t i=1;i<ang.size();i++) {
    switch (geom_.sphere()) {
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

    if (lonmin > longitude) lonmin = longitude;
    if (lonmax < longitude) lonmax = longitude;
    if (latmin > latitude) latmin = latitude;
    if (latmax < latitude) latmax = latitude;
  }

  geom_.setImageBounds(lonmin, lonmax, latmin, latmax);
}

void RenderArea::_findNewMaxResolution(Stomp::Map* stomp_map) {
  setMaxResolution(stomp_map->MaxResolution());
}


