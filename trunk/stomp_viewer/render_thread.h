// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the thread class for generating QPixmap from a
// Stomp::Map object.

#ifndef RENDER_THREAD_H
#define RENDER_THREAD_H

#include <QMutex>
#include <QThread>
#include <QWaitCondition>
#include <QPixmap>
#include <stomp.h>
#include "palette.h"
#include "render_geometry.h"

class Palette;
class RenderGeometry;
class RenderMapThread;
class RenderPointsThread;

class RenderMapThread : public QThread {
  Q_OBJECT

 public:
  RenderMapThread(QObject *parent = 0);
  ~RenderMapThread();

  void renderMap(Stomp::Map* stomp_map, RenderGeometry& geometry,
		 Palette& palette, uint32_t max_resolution,
		 bool antialiased, bool fill);

 signals:
  void renderedImage(const QImage& image);

 protected:
  void run();

 private:
  QMutex mutex;
  QWaitCondition condition;
  Stomp::Map* stomp_map_;
  RenderGeometry geom_;
  Palette palette_;
  uint32_t max_resolution_;
  bool antialiased_, fill_;
  bool restart;
  bool abort;
};

class RenderPointsThread : public QThread {
  Q_OBJECT

 public:
  RenderPointsThread(QObject *parent = 0);
  ~RenderPointsThread();

  void renderPoints(Stomp::WAngularVector* ang, RenderGeometry& geometry,
		    Palette& palette, bool antialiased);

 signals:
  void renderedImage(const QImage& image);

 protected:
  void run();

 private:
  QMutex mutex;
  QWaitCondition condition;
  Stomp::WAngularVector* ang_;
  RenderGeometry geom_;
  Palette palette_;
  bool antialiased_;
  bool restart;
  bool abort;
};

#endif

