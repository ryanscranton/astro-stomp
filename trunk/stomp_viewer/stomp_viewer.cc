#include <QtGui>
#include <QBoxLayout>

#include "render_area.h"
#include "stomp_viewer.h"

const int32_t IdRole = Qt::UserRole;

#define MESSAGE \
  Dialog::tr("<p>Message boxes have a caption, a text, "		\
	     "and any number of buttons, each with standard or custom texts." \
	     "<p>Click a button to close the message box. Pressing the Esc " \
	     "button will activate the detected escape button (if any).")

StompViewer::StompViewer(const QMap<QString, QSize> &customSizeHints,
                         QWidget *parent, Qt::WindowFlags flags) :
  QMainWindow(parent, flags) {
  renderArea = new RenderArea;
  renderArea->setImageBounds(-180.0, 180.0, -90.0, 90.0);
  renderArea->setWeightRange(0.0, 1.0);

  scrollArea = new QScrollArea;
  scrollArea->setBackgroundRole(QPalette::Dark);
  scrollArea->setWidgetResizable(false);
  scrollArea->setWidget(renderArea);
  errorMessageDialog = new QErrorMessage(this);

  dashLabel = new QLabel(tr(" - "));
  timesLabel = new QLabel(tr(" x "));

  pointsGroup = new QGroupBox(parent);
  // pointsGroup->setTitle(tr("Points Options"));
  pointsGroup->setMaximumWidth(400);

  setupMapIOGroup();
  setupMapCoordinatesGroup();
  setupMapAppearanceGroup();
  setupMapOtherGroup();

  setupPointsIOGroup();
  setupPointsAppearanceGroup();

  // Points group layout
  QVBoxLayout *pointsGroupLayout = new QVBoxLayout(pointsGroup);
  pointsGroupLayout->addWidget(pointsIOGroup);
  pointsGroupLayout->addWidget(pointsAppearanceGroup);
  pointsGroup->setLayout(pointsGroupLayout);

  // Our primary layout is a large rendering area on the left and a panel of
  // options for displaying the map on the right.
  setCentralWidget(scrollArea);

  coordinateSystemChanged();
  mapPaletteChanged();
  pointsPaletteChanged();
  resolutionChanged();
  mapWeightRangeChanged();
  pointsWeightRangeChanged();
  renderDimensionsChanged();

  setWindowTitle(tr("Stomp Map Viewer"));
  statusBar()->showMessage(tr("Ready"));
}

void StompViewer::setupMapIOGroup() {
  QToolBar* tool_bar = new QToolBar(this);
  tool_bar->setWindowTitle(tr("Map I/O"));
  addToolBar(tool_bar);

  QMenu *menu = new QMenu(tr("&File"), this);
  menuBar()->addMenu(menu);

  actionOpenMap = new QAction(QIcon(":/images/fileopen.png"),
			      tr("&Open New Map..."), this);
  actionOpenMap->setShortcut(QKeySequence::Open);
  connect(actionOpenMap, SIGNAL(triggered()),
	  this, SLOT(inputMapDialog()));
  tool_bar->addAction(actionOpenMap);
  menu->addAction(actionOpenMap);

  actionSaveImage = new QAction(QIcon(":/images/filesave.png"),
				tr("&Save to PNG..."), this);
  actionSaveImage->setShortcut(QKeySequence::Save);
  connect(actionSaveImage, SIGNAL(triggered()),
	  this, SLOT(saveImageDialog()));
  tool_bar->addAction(actionSaveImage);
  menu->addAction(actionSaveImage);

  menu->addSeparator();

  actionUpdateMap =
    new QAction(QIcon(":/images/editredo.png"), tr("&Update Map"), this);
  actionUpdateMap->setShortcut(Qt::CTRL + Qt::Key_U);
  connect(actionUpdateMap, SIGNAL(triggered()),
	  renderArea, SLOT(updatePixmap()));
  tool_bar->addAction(actionUpdateMap);
  menu->addAction(actionUpdateMap);

  actionAutoUpdate = new QAction(tr("Auto Upda&te"), this);
  actionAutoUpdate->setShortcut(Qt::CTRL + Qt::Key_T);
  actionAutoUpdate->setCheckable(true);
  connect(actionAutoUpdate, SIGNAL(toggled(bool)),
	  this, SLOT(autoUpdateToggled(bool)));
  menu->addAction(actionAutoUpdate);

  menu->addSeparator();

  actionLaunchPointsWindow = new QAction(QIcon(":/images/circle.png"),
					 tr("Overplot &points..."), this);
  actionLaunchPointsWindow->setShortcut(QKeySequence::Print);
  connect(actionLaunchPointsWindow, SIGNAL(triggered()),
	  this, SLOT(launchPointsWindow()));
  tool_bar->addAction(actionLaunchPointsWindow);
  menu->addAction(actionLaunchPointsWindow);

  // Connect any updates to the renderArea to a reset of the displayed map
  // parameters.
  connect(renderArea, SIGNAL(newMapParameters()),
	  this, SLOT(getNewMapParameters()));
}

void StompViewer::setupMapCoordinatesGroup() {
  QToolBar* tool_bar = new QToolBar(this);
  tool_bar->setWindowTitle(tr("Map Coordinates"));
  addToolBar(tool_bar);

  QMenu *menu = new QMenu(tr("&Coordinates"), this);
  menuBar()->addMenu(menu);

  actionAutoCoordinate = new QAction(tr("Auto Bounds"), this);
  actionAutoCoordinate->setShortcut(QKeySequence::Refresh);
  connect(actionAutoCoordinate, SIGNAL(triggered()),
	  this, SLOT(coordinateSystemChanged()));
  // tool_bar->addAction(actionOpenMap);
  menu->addAction(actionAutoCoordinate);
  menu->addSeparator();

  // Zoom buttons
  actionZoomIn =
    new QAction(QIcon(":/images/zoomin.png"), tr("Zoom &In"), this);
  actionZoomIn->setShortcut(QKeySequence::ZoomIn);
  connect(actionZoomIn, SIGNAL(triggered()), renderArea, SLOT(zoomIn()));
  connect(actionZoomIn, SIGNAL(triggered()),
	  this, SLOT(coordinateBoundsChanged()));
  tool_bar->addAction(actionZoomIn);
  menu->addAction(actionZoomIn);

  actionZoomOut =
    new QAction(QIcon(":/images/zoomout.png"), tr("Zoom &Out"), this);
  actionZoomOut->setShortcut(QKeySequence::ZoomOut);
  connect(actionZoomOut, SIGNAL(triggered()), renderArea, SLOT(zoomOut()));
  connect(actionZoomOut, SIGNAL(triggered()),
	  this, SLOT(coordinateBoundsChanged()));
  tool_bar->addAction(actionZoomOut);
  menu->addAction(actionZoomOut);

  // Begin system group -- within coordinates group
  systemGroup = new QGroupBox(tool_bar);
  systemGroup->setTitle(tr("Coordinate System"));

  coordinateComboBox = new QComboBox(systemGroup);
  coordinateComboBox->setWhatsThis(QString("Coordinate System"));
  coordinateComboBox->addItem(tr("Survey"), Stomp::AngularCoordinate::Survey);
  coordinateComboBox->addItem(tr("Equatorial"),
			      Stomp::AngularCoordinate::Equatorial);
  coordinateComboBox->addItem(tr("Galactic"),
			      Stomp::AngularCoordinate::Galactic);
  connect(coordinateComboBox, SIGNAL(activated(int)),
	  this, SLOT(coordinateSystemChanged()));

  QHBoxLayout *systemLayout = new QHBoxLayout(systemGroup);
  systemLayout->addWidget(coordinateComboBox);
  systemGroup->setLayout(systemLayout);
  tool_bar->addWidget(systemGroup);

  // Begin longitude group -- within coordinates group
  lonGroup = new QGroupBox(tool_bar);
  lonGroup->setTitle(tr("Longitude Range"));
  lonGroup->setMaximumWidth(200);

  lonMinLineEdit = new QLineEdit(lonGroup);
  lonMinLineEdit->setValidator(
    new QDoubleValidator(-180.0, 360.0, 6, lonMinLineEdit));
  lonMinLineEdit->setMaxLength(12);
  lonMaxLineEdit = new QLineEdit(lonGroup);
  lonMaxLineEdit->setValidator(
    new QDoubleValidator(-180.0, 360.0, 6, lonMaxLineEdit));
  lonMaxLineEdit->setMaxLength(12);
  connect(lonMinLineEdit, SIGNAL(editingFinished()),
	  this, SLOT(newLonMin()));
  connect(lonMaxLineEdit, SIGNAL(editingFinished()),
	  this, SLOT(newLonMax()));

  QHBoxLayout *lonLayout = new QHBoxLayout(lonGroup);
  lonLayout->addWidget(lonMinLineEdit);
  lonLayout->addWidget(dashLabel);
  lonLayout->addWidget(lonMaxLineEdit);
  lonGroup->setLayout(lonLayout);
  tool_bar->addWidget(lonGroup);

  // Begin longitude group -- within coordinates group
  latGroup = new QGroupBox(tool_bar);
  latGroup->setTitle(tr("Latitude Range"));
  latGroup->setMaximumWidth(200);

  latMinLineEdit = new QLineEdit(latGroup);
  latMinLineEdit->setValidator(
    new QDoubleValidator(-90.0, 90.0, 6, latMinLineEdit));
  latMinLineEdit->setMaxLength(10);
  latMaxLineEdit = new QLineEdit(latGroup);
  latMaxLineEdit->setValidator(
    new QDoubleValidator(-90.0, 90.0, 6, latMaxLineEdit));
  latMaxLineEdit->setMaxLength(10);
  connect(latMinLineEdit, SIGNAL(editingFinished()),
	  this, SLOT(newLatMin()));
  connect(latMaxLineEdit, SIGNAL(editingFinished()),
	  this, SLOT(newLatMax()));

  QHBoxLayout *latLayout = new QHBoxLayout(latGroup);
  latLayout->addWidget(latMinLineEdit);
  latLayout->addWidget(dashLabel);
  latLayout->addWidget(latMaxLineEdit);
  latGroup->setLayout(latLayout);
  tool_bar->addWidget(latGroup);

}

void StompViewer::setupMapAppearanceGroup() {
  QToolBar* tool_bar = new QToolBar(this);
  tool_bar->setWindowTitle(tr("Map Appearance"));
  tool_bar->setAllowedAreas(Qt::TopToolBarArea | Qt::BottomToolBarArea);
  addToolBarBreak(Qt::TopToolBarArea);
  addToolBar(tool_bar);

  // Begin resolution group -- within Appearance group
  QGroupBox *resolutionGroup = new QGroupBox(tool_bar);
  resolutionGroup->setTitle(tr("Maximum Pixel Resolution"));

  resolutionComboBox = new QComboBox(resolutionGroup);
  resolutionComboBox->addItem(tr("8"), 8);
  resolutionComboBox->addItem(tr("Default (2048)"), 2048);
  resolutionComboBox->addItem(tr("4"), 4);
  resolutionComboBox->addItem(tr("16"), 16);
  resolutionComboBox->addItem(tr("32"), 32);
  resolutionComboBox->addItem(tr("64"), 64);
  resolutionComboBox->addItem(tr("128"), 128);
  resolutionComboBox->addItem(tr("256"), 256);
  resolutionComboBox->addItem(tr("512"), 512);
  resolutionComboBox->addItem(tr("1024"), 1024);
  resolutionComboBox->addItem(tr("2048"), 2048);
  resolutionComboBox->addItem(tr("4096"), 4096);
  resolutionComboBox->addItem(tr("8192"), 8192);
  resolutionComboBox->addItem(tr("16384"), 16384);
  resolutionComboBox->addItem(tr("32768"), 32768);
  connect(resolutionComboBox, SIGNAL(activated(int)),
	  this, SLOT(resolutionChanged()));

  QHBoxLayout *resolutionLayout = new QHBoxLayout(resolutionGroup);
  resolutionLayout->addWidget(resolutionComboBox);
  resolutionGroup->setLayout(resolutionLayout);
  tool_bar->addWidget(resolutionGroup);

  // Begin palette group -- within Appearance group
  QGroupBox *paletteGroup = new QGroupBox(tool_bar);
  paletteGroup->setTitle(tr("Map Palette"));

  mapPaletteComboBox = new QComboBox(paletteGroup);
  mapPaletteComboBox->addItem(tr("Red Temperature"), Palette::RedTemperature);
  mapPaletteComboBox->addItem(tr("Green Temperature"),
			      Palette::GreenTemperature);
  mapPaletteComboBox->addItem(tr("Blue Temperature"), Palette::BlueTemperature);
  mapPaletteComboBox->addItem(tr("GrayScale"), Palette::GrayScale);
  mapPaletteComboBox->addItem(tr("Inverse GrayScale"),
			      Palette::InverseGrayScale);
  mapPaletteComboBox->addItem(tr("Rainbow"), Palette::Rainbow);
  mapPaletteComboBox->addItem(tr("Inverse Rainbow"), Palette::InverseRainbow);
  connect(mapPaletteComboBox, SIGNAL(activated(int)),
	  this, SLOT(mapPaletteChanged()));
  QHBoxLayout *paletteLayout = new QHBoxLayout(paletteGroup);
  paletteLayout->addWidget(mapPaletteComboBox);
  paletteGroup->setLayout(paletteLayout);
  tool_bar->addWidget(paletteGroup);

  // Begin weight group -- within Appearance group
  mapWeightGroup = new QGroupBox(tool_bar);
  mapWeightGroup->setTitle(tr("Weight Range"));
  mapWeightGroup->setMaximumWidth(200);

  mapWeightMinLineEdit = new QLineEdit(mapWeightGroup);
  mapWeightMinLineEdit->setValidator(
    new QDoubleValidator(-99999.0, 99999.0, 4, mapWeightMinLineEdit));
  mapWeightMinLineEdit->setMaxLength(8);
  mapWeightMaxLineEdit = new QLineEdit(mapWeightGroup);
  mapWeightMaxLineEdit->setValidator(
    new QDoubleValidator(-99999.0, 99999.0, 4, mapWeightMaxLineEdit));
  mapWeightMaxLineEdit->setMaxLength(8);
  connect(mapWeightMinLineEdit, SIGNAL(editingFinished()),
	  this, SLOT(newMapWeightMin()));
  connect(mapWeightMaxLineEdit, SIGNAL(editingFinished()),
	  this, SLOT(newMapWeightMax()));

  QHBoxLayout *weightLayout = new QHBoxLayout(mapWeightGroup);
  weightLayout->addWidget(mapWeightMinLineEdit);
  weightLayout->addWidget(dashLabel);
  weightLayout->addWidget(mapWeightMaxLineEdit);
  mapWeightGroup->setLayout(weightLayout);
  tool_bar->addWidget(mapWeightGroup);

  // Begin size group -- within Appearance group
  QGroupBox *sizeGroup = new QGroupBox(tool_bar);
  sizeGroup->setTitle(tr("Image Size (w x h)"));
  sizeGroup->setMaximumWidth(200);

  widthLineEdit = new QLineEdit(sizeGroup);
  widthLineEdit->setValidator(new QIntValidator(100, 32768, widthLineEdit));
  widthLineEdit->setMaxLength(5);
  heightLineEdit = new QLineEdit(sizeGroup);
  heightLineEdit->setValidator(new QIntValidator(100, 32768, heightLineEdit));
  heightLineEdit->setMaxLength(5);
  connect(widthLineEdit, SIGNAL(editingFinished()),
	  this, SLOT(newRenderWidth()));
  connect(heightLineEdit, SIGNAL(editingFinished()),
	  this, SLOT(newRenderHeight()));

  QHBoxLayout *sizeLayout = new QHBoxLayout(sizeGroup);
  sizeLayout->addWidget(widthLineEdit);
  sizeLayout->addWidget(timesLabel);
  sizeLayout->addWidget(heightLineEdit);
  sizeGroup->setLayout(sizeLayout);
  tool_bar->addWidget(sizeGroup);
}

void StompViewer::setupMapOtherGroup() {
  QMenu *menu = new QMenu(tr("&Options"), this);
  menuBar()->addMenu(menu);

  actionAntiAliasing = new QAction(tr("&Antialiasing"), this);
  actionAntiAliasing->setShortcut(Qt::CTRL + Qt::Key_A);
  actionAntiAliasing->setCheckable(true);
  connect(actionAntiAliasing, SIGNAL(toggled(bool)),
	  renderArea, SLOT(setAntialiased(bool)));
  menu->addAction(actionAntiAliasing);

  actionFillPixels = new QAction(tr("Fill Pixels"), this);
  actionFillPixels->setShortcut(Qt::CTRL + Qt::Key_D);
  actionFillPixels->setCheckable(true);
  connect(actionFillPixels, SIGNAL(toggled(bool)),
	  renderArea, SLOT(fillPolygons(bool)));
  menu->addAction(actionFillPixels);

  menu->addSeparator();

  actionAitoff = new QAction(tr("Aitoff-&Hammer Projection"), this);
  actionAitoff->setShortcut(Qt::CTRL + Qt::Key_H);
  actionAitoff->setCheckable(true);
  connect(actionAitoff, SIGNAL(toggled(bool)),
	  this, SLOT(aitoffToggled(bool)));
  menu->addAction(actionAitoff);

  actionFullSky = new QAction(tr("Fu&ll-Sky"), this);
  actionFullSky->setShortcut(Qt::CTRL + Qt::Key_L);
  actionFullSky->setCheckable(true);
  connect(actionFullSky, SIGNAL(toggled(bool)),
	  this, SLOT(fullSkyToggled(bool)));
  menu->addAction(actionFullSky);

  menu->addSeparator();

  actionShowCoordinates = new QAction(tr("Show Coor&dinates"), this);
  actionShowCoordinates->setShortcut(Qt::CTRL + Qt::Key_D);
  actionShowCoordinates->setCheckable(true);
  connect(actionShowCoordinates, SIGNAL(toggled(bool)),
	  this, SLOT(gridToggled()));
  menu->addAction(actionShowCoordinates);

  actionShowGrid = new QAction(tr("Show &Grid"), this);
  actionShowGrid->setShortcut(Qt::CTRL + Qt::Key_G);
  actionShowGrid->setCheckable(true);
  connect(actionShowGrid, SIGNAL(toggled(bool)),
	  this, SLOT(gridToggled()));
  menu->addAction(actionShowGrid);
}

void StompViewer::setupPointsIOGroup() {
  int32_t frameStyle = QFrame::Sunken | QFrame::Panel;

  // Now the section for inputing a new set of points.
  pointsIOGroup = new QGroupBox(pointsGroup);
  pointsIOGroup->setTitle(tr("Points I/O"));

  inputPointsButton = new QPushButton(tr("New Points..."), pointsIOGroup);
  currentPointsLabel = new QLabel(tr(""), pointsIOGroup);
  currentPointsLabel->setFrameStyle(frameStyle);

  pointsCoordinateComboBox = new QComboBox(pointsIOGroup);
  pointsCoordinateComboBox->addItem(tr("Survey Coordinates"),
				    Stomp::AngularCoordinate::Survey);
  pointsCoordinateComboBox->addItem(tr("Equatorial Coordinates"),
				    Stomp::AngularCoordinate::Equatorial);
  pointsCoordinateComboBox->addItem(tr("Galactic Coordinates"),
				    Stomp::AngularCoordinate::Galactic);

  weightedPointsCheckBox = new QCheckBox(tr("&Weighted Points"), pointsIOGroup);

  QGridLayout *ioLayout = new QGridLayout(pointsIOGroup);
  ioLayout->setColumnStretch(1, 1);
  ioLayout->setColumnMinimumWidth(1, 20);
  ioLayout->addWidget(inputPointsButton, 0, 0);
  ioLayout->addWidget(currentPointsLabel, 0, 1);
  ioLayout->addWidget(pointsCoordinateComboBox, 1, 0);
  ioLayout->addWidget(weightedPointsCheckBox, 1, 1);
  pointsIOGroup->setLayout(ioLayout);

  // Begin connections -- I/O group
  connect(inputPointsButton, SIGNAL(clicked()),
	  this, SLOT(inputPointsDialog()));
}

void StompViewer::setupPointsAppearanceGroup() {
  pointsAppearanceGroup = new QGroupBox(pointsGroup);
  pointsAppearanceGroup->setTitle(tr("Points Appearance"));

  showPointsCheckBox = new QCheckBox(tr("Show &Points"), pointsAppearanceGroup);
  filterPointsCheckBox = new QCheckBox(tr("Fil&ter Points"),
				       pointsAppearanceGroup);

  // Begin palette group -- within Appearance group
  pointsPaletteGroup = new QGroupBox(pointsAppearanceGroup);
  pointsPaletteGroup->setTitle(tr("Palette"));

  pointsPaletteComboBox = new QComboBox(pointsPaletteGroup);
  pointsPaletteComboBox->addItem(tr("Red Temperature"),
				 Palette::RedTemperature);
  pointsPaletteComboBox->addItem(tr("Green Temperature"),
				 Palette::GreenTemperature);
  pointsPaletteComboBox->addItem(tr("Blue Temperature"),
				 Palette::BlueTemperature);
  pointsPaletteComboBox->addItem(tr("GrayScale"), Palette::GrayScale);
  pointsPaletteComboBox->addItem(tr("Inverse GrayScale"),
				 Palette::InverseGrayScale);
  pointsPaletteComboBox->addItem(tr("Rainbow"), Palette::Rainbow);
  pointsPaletteComboBox->addItem(tr("Inverse Rainbow"),
				 Palette::InverseRainbow);

  QHBoxLayout *paletteLayout = new QHBoxLayout(pointsPaletteGroup);
  paletteLayout->addWidget(pointsPaletteComboBox);
  pointsPaletteGroup->setLayout(paletteLayout);

  // Begin weight group -- within Appearance group
  pointsWeightGroup = new QGroupBox(pointsAppearanceGroup);
  pointsWeightGroup->setTitle(tr("Weight Range"));

  pointsWeightMinLineEdit = new QLineEdit(pointsWeightGroup);
  pointsWeightMinLineEdit->setValidator(
    new QDoubleValidator(-99999.0, 99999.0, 4, pointsWeightMinLineEdit));
  pointsWeightMinLineEdit->setMaxLength(8);
  pointsWeightMaxLineEdit = new QLineEdit(pointsWeightGroup);
  pointsWeightMaxLineEdit->setValidator(
    new QDoubleValidator(-99999.0, 99999.0, 4, pointsWeightMaxLineEdit));
  pointsWeightMaxLineEdit->setMaxLength(8);

  QHBoxLayout *weightLayout = new QHBoxLayout(pointsWeightGroup);
  weightLayout->addWidget(pointsWeightMinLineEdit);
  weightLayout->addWidget(dashLabel);
  weightLayout->addWidget(pointsWeightMaxLineEdit);
  pointsWeightGroup->setLayout(weightLayout);

  QGridLayout *appearanceLayout = new QGridLayout(pointsAppearanceGroup);
  appearanceLayout->addWidget(showPointsCheckBox, 0, 0);
  appearanceLayout->addWidget(filterPointsCheckBox, 0, 1);
  appearanceLayout->addWidget(pointsWeightGroup, 1, 0);
  appearanceLayout->addWidget(pointsPaletteGroup, 1, 1);
  pointsAppearanceGroup->setLayout(appearanceLayout);

  // Appearance group connections
  connect(showPointsCheckBox, SIGNAL(toggled(bool)),
	  this, SLOT(showPointsToggled()));
  connect(filterPointsCheckBox, SIGNAL(toggled(bool)),
	  this, SLOT(filterPointsToggled()));
  connect(pointsWeightMinLineEdit, SIGNAL(editingFinished()),
	  this, SLOT(newPointsWeightMin()));
  connect(pointsWeightMaxLineEdit, SIGNAL(editingFinished()),
	  this, SLOT(newPointsWeightMax()));
  connect(pointsPaletteComboBox, SIGNAL(activated(int)),
	  this, SLOT(pointsPaletteChanged()));
}

void StompViewer::inputMapDialog() {
  // Before we begin, we disable updates on the renderArea so that the map
  // doesn't re-draw due to our dialog bog
  renderArea->setUpdatesEnabled(false);

  // launch our dialog box
  QFileDialog::Options options = 0;
#ifdef Q_WS_MAC
  options = QFileDialog::ReadOnly;
#endif
  QString selectedFilter;
  QString fileName =
    QFileDialog::getOpenFileName(this,
				 tr("Open a new Stomp Map"),
                                 QDir::currentPath(),
                                 tr("All Files (*)"),
                                 &selectedFilter,
                                 options);
  if (!fileName.isEmpty()) {
    statusBar()->showMessage(QString("Reading %1").arg(fileName));
    if (renderArea->readNewMap(fileName)) {
      if (fileName.contains("/")) {
	statusBar()->showMessage(QString("Loaded %1").
				 arg(fileName.section("/", -1)));
	setWindowTitle(tr("Stomp Map Viewer: ")+
		       QString(QString("%1").
			       arg(fileName.section("/", -1))));
      } else {
	statusBar()->showMessage(QString("Loaded %1").arg(fileName));
	setWindowTitle(tr("Stomp Map Viewer: ")+
		       QString(QString("%1").arg(fileName)));
      }
    }
  }

  // And re-enable updates on renderArea.
  renderArea->setUpdatesEnabled(true);
}

void StompViewer::saveImageDialog() {
  // Before we begin, we disable updates on the renderArea so that the map
  // doesn't re-draw due to our dialog bog
  renderArea->setUpdatesEnabled(false);

  // launch our dialog box
  QFileDialog::Options options;
  QString selectedFilter;
  QString fileName =
    QFileDialog::getSaveFileName(this,
				 tr("Save image as a PNG"),
                                 QDir::currentPath(),
                                 tr("Image files (*.png, *.PNG)"),
                                 &selectedFilter,
                                 options);
  if (!fileName.isEmpty()) {
    statusBar()->showMessage(QString("Writing to %1").arg(fileName));
    if (renderArea->writeToPng(fileName)) {
      if (fileName.contains("/")) {
	statusBar()->showMessage(QString("Wrote image to  %1").
				 arg(fileName.section("/", -1)));
      } else {
	statusBar()->showMessage(QString("Wrote image to %1").arg(fileName));
      }
    }
  }

  // And re-enable updates on renderArea.
  renderArea->setUpdatesEnabled(true);
}

void StompViewer::coordinateBoundsChanged() {
  lonMinLineEdit->setText(QString("%1").arg(renderArea->longitudeMin()));
  lonMaxLineEdit->setText(QString("%1").arg(renderArea->longitudeMax()));
  latMinLineEdit->setText(QString("%1").arg(renderArea->latitudeMin()));
  latMaxLineEdit->setText(QString("%1").arg(renderArea->latitudeMax()));
}

void StompViewer::coordinateSystemChanged() {
  Stomp::AngularCoordinate::Sphere sphere =
    Stomp::AngularCoordinate::Sphere(coordinateComboBox->itemData(
				       coordinateComboBox->currentIndex(),
				       IdRole).toInt());
  renderArea->setCoordinateSystem(sphere);
  coordinateBoundsChanged();
}

void StompViewer::newLonMin() {
  double lonmin = lonMinLineEdit->text().toDouble();
  renderArea->setImageBounds(lonmin,
			     renderArea->longitudeMax(),
			     renderArea->latitudeMin(),
			     renderArea->latitudeMax());
}

void StompViewer::newLonMax() {
  double lonmax = lonMaxLineEdit->text().toDouble();
  renderArea->setImageBounds(renderArea->longitudeMin(),
			     lonmax,
			     renderArea->latitudeMin(),
			     renderArea->latitudeMax());
}

void StompViewer::newLatMin() {
  double latmin = latMinLineEdit->text().toDouble();
  renderArea->setImageBounds(renderArea->longitudeMin(),
			     renderArea->longitudeMax(),
			     latmin,
			     renderArea->latitudeMax());
}

void StompViewer::newLatMax() {
  double latmax = latMaxLineEdit->text().toDouble();
  renderArea->setImageBounds(renderArea->longitudeMin(),
			     renderArea->longitudeMax(),
			     renderArea->latitudeMin(),
			     latmax);
}

void StompViewer::mapPaletteChanged() {
  Palette::PaletteType palette_type =
    Palette::PaletteType(mapPaletteComboBox->itemData(
			   mapPaletteComboBox->currentIndex(),
			   IdRole).toInt());
  renderArea->setPalette(palette_type);
}

void StompViewer::resolutionChanged() {
  int resolution =
    resolutionComboBox->itemData(resolutionComboBox->currentIndex(),
				 IdRole).toInt();
  renderArea->setMaxResolution(static_cast<uint16_t>(resolution));
}

void StompViewer::newMapWeightMin() {
  double weight_min = mapWeightMinLineEdit->text().toDouble();
  double weight_max = renderArea->weightMax();

  if (weight_min > weight_max) {
    double tmp_weight = weight_min;

    weight_min = weight_max;
    weight_max = tmp_weight;
  }

  renderArea->setWeightRange(weight_min, weight_max);
  mapWeightRangeChanged();
}

void StompViewer::newMapWeightMax() {
  double weight_max = mapWeightMaxLineEdit->text().toDouble();
  double weight_min = renderArea->weightMin();

  if (weight_min > weight_max) {
    double tmp_weight = weight_min;

    weight_min = weight_max;
    weight_max = tmp_weight;
  }

  renderArea->setWeightRange(weight_min, weight_max);
  mapWeightRangeChanged();
}

void StompViewer::mapWeightRangeChanged() {
  mapWeightMinLineEdit->setText(QString("%1").arg(renderArea->weightMin()));
  mapWeightMaxLineEdit->setText(QString("%1").arg(renderArea->weightMax()));
}

void StompViewer::pointsWeightRangeChanged() {
  pointsWeightMinLineEdit->setText(QString("%1").
				   arg(renderArea->pointsWeightMin()));
  pointsWeightMaxLineEdit->setText(QString("%1").
				   arg(renderArea->pointsWeightMax()));
}

void StompViewer::newRenderWidth() {
  int render_width = widthLineEdit->text().toInt();
  int render_height = renderArea->height();

  if (render_width != renderArea->width()) {
    // int panel_width = width() - renderArea->width();

    renderArea->setFixedSize(render_width, render_height);
    // resize(render_width+panel_width, height());
    renderDimensionsChanged();
  }
}

void StompViewer::newRenderHeight() {
  int render_width = renderArea->width();
  int render_height = heightLineEdit->text().toInt();

  if (render_height != renderArea->height()) {
    renderArea->setFixedSize(render_width, render_height);
    // if (render_height > height())
    // resize(width(), render_height);
    renderDimensionsChanged();
  }
}

void StompViewer::renderDimensionsChanged() {
  widthLineEdit->setText(QString("%1").arg(renderArea->width()));
  heightLineEdit->setText(QString("%1").arg(renderArea->height()));
}

void StompViewer::fullSkyToggled(bool full_sky) {
  renderArea->setFullSky(full_sky);
  coordinateBoundsChanged();

  // if we're in full-sky mode, then we need to disable the buttons for
  // changing the coordinate bounds.
  if (full_sky) {
    actionAutoCoordinate->blockSignals(true);
    lonGroup->setDisabled(true);
    latGroup->setDisabled(true);
  } else {
    // Likewise, if we've disabled full-sky mode, then we need to re-activate
    // the buttons.
    actionAutoCoordinate->blockSignals(false);
    lonGroup->setDisabled(false);
    latGroup->setDisabled(false);

    // If we're no longer in full-sky mode, then we need to uncheck the Aitoff
    // box if it's checked.
    if (actionAitoff->isChecked()) actionAitoff->setChecked(false);
  }
}

void StompViewer::aitoffToggled(bool use_aitoff) {
  renderArea->useAitoffProjection(use_aitoff);

  // If we've turned on the Aitoff projection, then we need to also set our
  // coordinate bounds to full-sky.
  if (use_aitoff) actionFullSky->setChecked(true);
}

void StompViewer::gridToggled() {
  renderArea->setDisplayCoordinates(actionShowCoordinates->isChecked());
  renderArea->setDisplayGrid(actionShowGrid->isChecked());

  renderArea->updateGrid();
}

void StompViewer::autoUpdateToggled(bool auto_update) {
  renderArea->setAutoUpdate(auto_update);

  // Need to disable the update button if we've gone to auto-updating.  Or
  // vice-versa.
  if (auto_update) {
    actionUpdateMap->blockSignals(true);
  } else {
    actionUpdateMap->blockSignals(false);
  }
}

void StompViewer::getNewMapParameters() {
  bool set_to_auto_update = false;
  if (renderArea->autoUpdating()) {
    set_to_auto_update = true;
    renderArea->setAutoUpdate(false);
  }

  getNewWeightRange();
  getNewMaxResolution();
  coordinateBoundsChanged();
  if (renderArea->fullSky()) fullSkyToggled(true);

  renderArea->setAutoUpdate(set_to_auto_update);
  if (renderArea->autoUpdating()) renderArea->update();
}

void StompViewer::getNewWeightRange() {
  double weight_min = renderArea->weightMin();
  double weight_max = renderArea->weightMax();

  if (Stomp::DoubleEQ(weight_min, weight_max)) {
    renderArea->setWeightRange(0.0, weight_max);
  }
  mapWeightRangeChanged();

  weight_min = renderArea->pointsWeightMin();
  weight_max = renderArea->pointsWeightMax();

  if (Stomp::DoubleEQ(weight_min, weight_max)) {
    renderArea->setPointsWeightRange(0.0, weight_max);
  }
  pointsWeightRangeChanged();
}

void StompViewer::getNewMaxResolution() {
  int resolution = renderArea->maxResolution();

  // Iterate over the possible resolution values until we find one that
  // matches the current maximum resolution.
  int idx = 0;
  resolutionComboBox->setCurrentIndex(idx);

  while (resolution != resolutionComboBox->itemData(
	   resolutionComboBox->currentIndex(), IdRole).toInt() &&
	 idx < resolutionComboBox->count()) {
    idx++;
    resolutionComboBox->setCurrentIndex(idx);
  }
}

void StompViewer::launchPointsWindow() {
  QDockWidget *pointsDockWidget = new QDockWidget(tr("Overplot Points"), this);
  pointsDockWidget->setFeatures(QDockWidget::AllDockWidgetFeatures);
  pointsDockWidget->setAllowedAreas(Qt::LeftDockWidgetArea);
  pointsDockWidget->setWidget(pointsGroup);
  addDockWidget(Qt::LeftDockWidgetArea, pointsDockWidget);
  pointsDockWidget->setFloating(true);
}

void StompViewer::inputPointsDialog() {
  // Before we begin, we disable updates on the renderArea so that the map
  // doesn't re-draw due to our dialog bog
  renderArea->setUpdatesEnabled(false);

  // Grab our input parameters from the I/O elements
  Stomp::AngularCoordinate::Sphere sphere =
    Stomp::AngularCoordinate::Sphere(pointsCoordinateComboBox->itemData(
				       pointsCoordinateComboBox->currentIndex(),
				       IdRole).toInt());
  bool weighted_points = weightedPointsCheckBox->isChecked();

  // launch our dialog box
  QFileDialog::Options options = 0;
#ifdef Q_WS_MAC
  options = QFileDialog::ReadOnly;
#endif
  QString selectedFilter;
  QString fileName =
    QFileDialog::getOpenFileName(this,
				 tr("Open a file containing input points"),
                                 QDir::currentPath(),
                                 tr("All Files (*)"),
                                 &selectedFilter,
                                 options);
  if (!fileName.isEmpty()) {
    renderArea->clearPoints();
    if (!showPointsCheckBox->isChecked()) showPointsCheckBox->setChecked(true);
    statusBar()->showMessage(QString("Reading %1").arg(fileName));
    if (renderArea->readNewPointsFile(fileName, sphere, weighted_points)) {
      if (fileName.contains("/")) {
	currentPointsLabel->setText(fileName.section("/", -1));
      } else {
	currentPointsLabel->setText(fileName);
      }
    }
  }

  // And re-enable updates on renderArea.
  renderArea->setUpdatesEnabled(true);
}

void StompViewer::pointsPaletteChanged() {
  Palette::PaletteType palette_type =
    Palette::PaletteType(pointsPaletteComboBox->itemData(
			   pointsPaletteComboBox->currentIndex(),
			   IdRole).toInt());
  renderArea->setPointsPalette(palette_type);
}

void StompViewer::newPointsWeightMin() {
  double weight_min = pointsWeightMinLineEdit->text().toDouble();
  double weight_max = renderArea->pointsWeightMax();

  if (weight_min > weight_max) {
    double tmp_weight = weight_min;

    weight_min = weight_max;
    weight_max = tmp_weight;
  }

  renderArea->setPointsWeightRange(weight_min, weight_max);
  pointsWeightRangeChanged();
}

void StompViewer::newPointsWeightMax() {
  double weight_max = pointsWeightMaxLineEdit->text().toDouble();
  double weight_min = renderArea->pointsWeightMin();

  if (weight_min > weight_max) {
    double tmp_weight = weight_min;

    weight_min = weight_max;
    weight_max = tmp_weight;
  }

  renderArea->setPointsWeightRange(weight_min, weight_max);
  pointsWeightRangeChanged();
}

void StompViewer::showPointsToggled() {
  renderArea->setDisplayPoints(showPointsCheckBox->isChecked());
  renderArea->setFilterPoints(filterPointsCheckBox->isChecked());

  renderArea->updatePoints();
}

void StompViewer::filterPointsToggled() {
  renderArea->setDisplayPoints(showPointsCheckBox->isChecked());
  renderArea->setFilterPoints(filterPointsCheckBox->isChecked());

  renderArea->updatePoints();
}
