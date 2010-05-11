#include <QtGui>
#include "render_geometry.h"

RenderGeometry::RenderGeometry() {
  buffer_top_ = buffer_bottom_ = buffer_left_ = buffer_right_ = 0;
  aitoff_projection_ = false;

  width_ = 0;
  height_ = 0;

  setImageBounds(0.0, 360.0, -90.0, 90.0,
		 Stomp::AngularCoordinate::Equatorial);
}

RenderGeometry::RenderGeometry(double longitude_min, double longitude_max,
			       double latitude_min, double latitude_max,
			       Stomp::AngularCoordinate::Sphere sphere,
			       int width, int height) {
  buffer_top_ = buffer_bottom_ = buffer_left_ = buffer_right_ = 0;
  aitoff_projection_ = false;

  width_ = width;
  height_ = height;

  setImageBounds(longitude_min, longitude_max, latitude_min, latitude_max,
		 sphere);
}

double RenderGeometry::longitudeMin() {
  return lonmin_;
}

double RenderGeometry::longitudeMax() {
  return lonmax_;
}

double RenderGeometry::latitudeMin() {
  return latmin_;
}

double RenderGeometry::latitudeMax() {
  return latmax_;
}

Stomp::AngularCoordinate::Sphere RenderGeometry::sphere() {
  return sphere_;
}

int RenderGeometry::height() {
  return height_;
}

int RenderGeometry::width() {
  return width_;
}

bool RenderGeometry::aitoffProjection() {
  return aitoff_projection_;
}

bool RenderGeometry::fullSky() {
  return full_sky_;
}

void RenderGeometry::setVerticalBuffer(uint16_t buffer_pixels) {
  buffer_top_ = buffer_pixels;
  buffer_bottom_ = buffer_pixels;
}

void RenderGeometry::setHorizontalBuffer(uint16_t buffer_pixels) {
  buffer_left_ = buffer_pixels;
  buffer_right_ = buffer_pixels;
}

void RenderGeometry::setTopBuffer(uint16_t buffer_pixels) {
  buffer_top_ = buffer_pixels;
}

void RenderGeometry::setBottomBuffer(uint16_t buffer_pixels) {
  buffer_bottom_ = buffer_pixels;
}

void RenderGeometry::setLeftBuffer(uint16_t buffer_pixels) {
  buffer_left_ = buffer_pixels;
}

void RenderGeometry::setRightBuffer(uint16_t buffer_pixels) {
  buffer_right_ = buffer_pixels;
}

void RenderGeometry::setImageBounds(double lonmin, double lonmax,
				    double latmin, double latmax) {
  latmin_ = latmin;
  latmax_ = latmax;
  lonmin_ = lonmin;
  lonmax_ = lonmax;
  checkBounds();
}

void RenderGeometry::setImageBounds(double lonmin, double lonmax,
				    double latmin, double latmax,
				    Stomp::AngularCoordinate::Sphere sphere) {
  latmin_ = latmin;
  latmax_ = latmax;
  lonmin_ = lonmin;
  lonmax_ = lonmax;
  sphere_ = sphere;
  checkBounds();
}

void RenderGeometry::setCoordinateSystem(Stomp::AngularCoordinate::Sphere sphere) {
  sphere_ = sphere;
}

void RenderGeometry::setHeight(int height) {
  height_ = height;
}

void RenderGeometry::setWidth(int width) {
  width_ = width;
}

void RenderGeometry::setSize(const QSize& new_size) {
  height_ = new_size.height();
  width_ = new_size.width();
}

void RenderGeometry::useAitoffProjection(bool use_aitoff_projection) {
  aitoff_projection_ = use_aitoff_projection;
}

void RenderGeometry::setFullSky(bool full_sky) {
  if (full_sky) {
    double latmin = -90.0, latmax = 90.0;
    double lonmin = 0.0, lonmax = 360.0;
    if (sphere() == Stomp::AngularCoordinate::Survey) {
      lonmin = -180.0;
      lonmax = 180.0;
    }

    setImageBounds(lonmin, lonmax, latmin, latmax);
    full_sky_ = true;
  } else {
    full_sky_ = false;
  }
}

void RenderGeometry::checkBounds() {
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

void RenderGeometry::enforceBounds(double& lon, double& lat) {
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

void RenderGeometry::pixelToPoint(uint32_t pixel_x, uint32_t pixel_y,
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

  enforceBounds(lon_center, lat_center);

  if (aitoff_projection_) cartesianToAitoff(lon_center, lat_center);

  point.setX(lonToX(lon_center));
  point.setY(latToY(lat_center));
}

void RenderGeometry::pixelToPolygon(uint32_t pixel_x, uint32_t pixel_y,
				    uint32_t resolution,
				    uint32_t max_resolution,
				    QPolygonF& polygon) {
  if (!polygon.empty()) polygon.clear();

  uint32_t ratio = max_resolution/resolution;

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
    enforceBounds(lon, lat);

    if (aitoff_projection_) cartesianToAitoff(lon, lat);
    polygon << QPointF(lonToX(lon), latToY(lat));
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
    enforceBounds(lon, lat);

    if (aitoff_projection_) cartesianToAitoff(lon, lat);
    polygon << QPointF(lonToX(lon), latToY(lat));
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
    enforceBounds(lon, lat);

    if (aitoff_projection_) cartesianToAitoff(lon, lat);
    polygon << QPointF(lonToX(lon), latToY(lat));
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
    enforceBounds(lon, lat);

    if (aitoff_projection_) cartesianToAitoff(lon, lat);
    polygon << QPointF(lonToX(lon), latToY(lat));
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
  enforceBounds(lon, lat);

  if (aitoff_projection_) cartesianToAitoff(lon, lat);
  polygon << QPointF(lonToX(lon), latToY(lat));
}

void RenderGeometry::splitPixelToPolygons(uint32_t pixel_x, uint32_t pixel_y,
					  uint32_t resolution,
					  uint32_t max_resolution,
					  QPolygonF& left_polygon,
					  QPolygonF& right_polygon) {
  if (!left_polygon.empty()) left_polygon.clear();
  if (!right_polygon.empty()) right_polygon.clear();

  uint32_t ratio = max_resolution/resolution;

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
    enforceBounds(left_lon, lat);
    enforceBounds(right_lon, lat);
    if (aitoff_projection_) {
      cartesianToAitoff(left_lon, lat);
      cartesianToAitoff(right_lon, lat);
    }

    left_polygon << QPointF(lonToX(left_lon), latToY(lat));
    right_polygon << QPointF(lonToX(right_lon), latToY(lat));
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
    enforceBounds(left_lon, lat);
    enforceBounds(right_lon, lat);
    if (aitoff_projection_) {
      cartesianToAitoff(left_lon, lat);
      cartesianToAitoff(right_lon, lat);
    }

    left_polygon << QPointF(lonToX(left_lon), latToY(lat));
    right_polygon << QPointF(lonToX(right_lon), latToY(lat));
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
    enforceBounds(left_lon, lat);
    enforceBounds(right_lon, lat);
    if (aitoff_projection_) {
      cartesianToAitoff(left_lon, lat);
      cartesianToAitoff(right_lon, lat);
    }

    left_polygon << QPointF(lonToX(left_lon), latToY(lat));
    right_polygon << QPointF(lonToX(right_lon), latToY(lat));
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
    enforceBounds(left_lon, lat);
    enforceBounds(right_lon, lat);
    if (aitoff_projection_) {
      cartesianToAitoff(left_lon, lat);
      cartesianToAitoff(right_lon, lat);
    }

    left_polygon << QPointF(lonToX(left_lon), latToY(lat));
    right_polygon << QPointF(lonToX(right_lon), latToY(lat));
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
  enforceBounds(left_lon, lat);
  enforceBounds(right_lon, lat);
  if (aitoff_projection_) {
    cartesianToAitoff(left_lon, lat);
    cartesianToAitoff(right_lon, lat);
  }

  left_polygon << QPointF(lonToX(left_lon), latToY(lat));
  right_polygon << QPointF(lonToX(right_lon), latToY(lat));
}

bool RenderGeometry::angToPoint(Stomp::WeightedAngularCoordinate& ang,
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

    if (aitoff_projection_) cartesianToAitoff(longitude, latitude);

    point.setX(lonToX(longitude));
    point.setY(latToY(latitude));
  }

  return inside_bounds;
}

qreal RenderGeometry::lonToX(double longitude) {
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

  double pixel_x = range*n_pixel/lonRange();
  if (Stomp::DoubleLT(pixel_x, 0.0)) pixel_x = 0.0;
  if (Stomp::DoubleGE(pixel_x, n_pixel))
    pixel_x = n_pixel;

  return static_cast<qreal>(pixel_x + buffer_left_);
}

double RenderGeometry::xToLon(qreal pixel_x) {
  double n_pixel = static_cast<double>(width() - buffer_left_ - buffer_right_);
  pixel_x -= static_cast<qreal>(buffer_left_);

  double range = static_cast<double>(lonRange()*(pixel_x/n_pixel));

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

qreal RenderGeometry::latToY(double latitude) {
  double n_pixel = static_cast<double>(height() - buffer_top_ - buffer_bottom_);
  double range = latmax_ - latitude;

  double pixel_y = range*n_pixel/latRange();
  if (Stomp::DoubleLT(pixel_y, 0.0)) pixel_y = 0.0;
  if (Stomp::DoubleGE(pixel_y, n_pixel)) pixel_y = n_pixel;

  return static_cast<qreal>(pixel_y + buffer_top_);
}

double RenderGeometry::yToLat(qreal pixel_y) {
  double n_pixel = static_cast<double>(height() - buffer_top_ - buffer_bottom_);
  pixel_y -= static_cast<qreal>(buffer_top_);

  double range = static_cast<double>(latRange()*(pixel_y/n_pixel));
  double latitude = latmax_ - range;

  if (Stomp::DoubleLE(latitude, latmin_)) latitude = latmin_;
  if (Stomp::DoubleGE(latitude, latmax_)) latitude = latmax_;

  return latitude;
}

double RenderGeometry::lonRange() {
  return (continuous_longitude_ ? lonmax_ - lonmin_ :
	  360.0 - lonmin_ + lonmax_);
}

double RenderGeometry::latRange() {
  return latmax_ - latmin_;
}

double RenderGeometry::lonCenter() {
  double lon_center = 0.5*(lonmin_ + lonmax_);
  if (!continuous_longitude_) {
    lon_center = lonmax_ + 0.5*lonRange();
    if (sphere_ == Stomp::AngularCoordinate::Survey) {
      if (Stomp::DoubleGE(lon_center, 180.0)) lon_center -= 360.0;
    } else {
      if (Stomp::DoubleGE(lon_center, 360.0)) lon_center -= 360.0;
    }
  }

  return lon_center;
}

double RenderGeometry::latCenter() {
  return 0.5*(latmin_ + latmax_);
}

double RenderGeometry::renderPixelArea() {
  double n_pixel =
    static_cast<double>((width() - buffer_left_ - buffer_right_)*
			(height() - buffer_top_ - buffer_bottom_));

  double cartesian_area = lonRange()*latRange();

  return cartesian_area/n_pixel;
}

void RenderGeometry::cartesianToAitoff(double& longitude, double& latitude) {
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

