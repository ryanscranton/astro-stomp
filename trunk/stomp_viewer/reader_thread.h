// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the thread class for reading a Stomp::Map file and
// generating softened versions of it for rendering.

#ifndef READER_THREAD_H
#define READER_THREAD_H

#include <QMutex>
#include <QThread>
#include <QWaitCondition>
#include <stomp.h>

class ReadMapThread;

class ReadMapThread : public QThread {
  Q_OBJECT

 public:
  ReadMapThread(QObject *parent = 0);
  ~ReadMapThread();

  void readMap(const std::string& map_file_name);
  void clearMaps();

 signals:
  void newBaseMap(Stomp::Map* base_map);
  void newSoftenedMap(uint8_t level, Stomp::Map* softened_map);
  void readerProgress(int progress);

 protected:
  void run();

 private:
  QMutex mutex;
  QWaitCondition condition;
  std::string map_file_name_;
  Stomp::Map* base_map_;
  std::map<uint8_t, Stomp::Map*> softened_maps_;
  bool good_base_map_;
  uint8_t softened_level_;
  bool restart;
  bool abort;
};

#endif

