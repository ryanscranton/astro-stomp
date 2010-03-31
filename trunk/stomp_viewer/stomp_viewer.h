#ifndef STOMP_VIEWER_H
#define STOMP_VIEWER_H

#include <QMainWindow>
#include <QWidget>
#include <QAction>
#include <QToolBar>
#include <QMenu>
#include <QMenuBar>
#include <QDialog>
#include <QFrame>
#include <QGroupBox>
#include <QLineEdit>
#include <QPushButton>

class QCheckBox;
class QComboBox;
class QScrollArea;
class QErrorMessage;
class QGroupBox;
class QLabel;
class QLineEdit;
class QPushButton;
class RenderArea;

class StompViewer : public QMainWindow {
  Q_OBJECT

 public:
  StompViewer(const QMap<QString, QSize> &customSizeHints,
	      QWidget *parent = 0, Qt::WindowFlags flags = 0);

 public slots:
  // MAP I/O SLOTS
  //
  // First, the slot to launch a dialog box to read in a new map.  Provided this
  // works, the renderArea will send back a signal alerting the StompViewer
  // to update its display information.
  void inputMapDialog();

  // Second, we've got the slot that handles requests to write the renderArea
  // out to a png.  This will launch a dialog to get the output file name
  // and send the appropriate commands to renderArea.
  void saveImageDialog();

  // MAP COORDINATE SLOTS
  //
  // Any time that the coordinates for the map might have been changed,
  // we need to update the values displayed in the right side panel.
  void coordinateBoundsChanged();

  // This handles the QComboBox that specifies the coordinate system that we're
  // showing at the moment.
  void coordinateSystemChanged();

  //  These four methods are tied to the QPushButtons that allow for changing
  // the lat-lon bounds of the area displayed.
  void newLonMin();
  void newLonMax();
  void newLatMin();
  void newLatMax();

  // MAP APPEARANCE SLOTS
  //
  // First, the slots that handles changes to the weight bounds that we're
  // using.  These operate just like the coordinate bounds slots.
  void newMapWeightMin();
  void newMapWeightMax();
  void mapWeightRangeChanged();

  // Now, some methods to process changes to the renderArea dimensions.
  void newRenderWidth();
  void newRenderHeight();
  void renderDimensionsChanged();

  // Next, a slot to process changes to the image palette.
  void mapPaletteChanged();

  // And finally, a slot to process the QComboBox controlling the maximum
  // pixel resolution to use when generating the map image.
  void resolutionChanged();

  //  CHECKBOX SLOTS
  //
  // Slot to process requests that the map be switched to full-sky.
  // When this is checked, the buttons for manually adjusting the map bounds
  // are disabled and vice versa.
  void fullSkyToggled(bool full_sky);

  // Slot to process toggling of the Aitoff-Hammer checkbox.  If this is checked
  // then the full-sky box needs to be checked as well.
  void aitoffToggled(bool use_aitoff);

  // Slot to process requests that the coordinates or coordinate grid be
  // turned on or off.  Since the mechanics are similar in both cases, we'll
  // handle both instances with a single method.
  void gridToggled();

  //  MAP UPDATE SLOTS
  //
  // Just one new slot here, the one that processes a request to auto-update
  // the image in response to changes in the right panel.  If this is checked
  // then the update button is disabled.
  void autoUpdateToggled(bool auto_update);

  // NEW MAP SLOTS
  void getNewMapParameters();
  void getNewWeightRange();
  void getNewMaxResolution();

  // POINTS I/O SLOTS
  //
  // The slot to open the dialog for inputing a new set of point coordinates.
  void launchPointsWindow();
  void inputPointsDialog();

  // POINTS APPEARANCE SLOTS
  //
  // First, the slots that handles changes to the weight bounds that we're
  // using.  These operate just like the coordinate bounds slots.
  void newPointsWeightMin();
  void newPointsWeightMax();
  void pointsWeightRangeChanged();

  // Next, a slot to process changes to the image palette.
  void pointsPaletteChanged();

  // And two slots to handle the check boxes to show and/or filter the points
  // against the map.
  void showPointsToggled();
  void filterPointsToggled();

 private:
  void setupMapIOGroup();
  void setupMapCoordinatesGroup();
  void setupMapAppearanceGroup();
  void setupMapOtherGroup();
  void setupPointsIOGroup();
  void setupPointsAppearanceGroup();

  RenderArea *renderArea;
  QScrollArea *scrollArea;

  QGroupBox *pointsGroup;
  QAction *actionLaunchPointsWindow;

  QAction *actionOpenMap;
  QAction *actionSaveImage;
  QAction *actionUpdateMap;
  QAction *actionAutoUpdate;

  QGroupBox *systemGroup;
  QComboBox *coordinateComboBox;
  QGroupBox *lonGroup;
  QGroupBox *latGroup;
  QLineEdit *lonMinLineEdit;
  QLineEdit *lonMaxLineEdit;
  QLineEdit *latMinLineEdit;
  QLineEdit *latMaxLineEdit;
  QAction *actionAutoCoordinate;
  QAction *actionZoomIn;
  QAction *actionZoomOut;

  QGroupBox *mapWeightGroup;
  QLineEdit *mapWeightMinLineEdit;
  QLineEdit *mapWeightMaxLineEdit;
  QGroupBox *sizeGroup;
  QLineEdit *widthLineEdit;
  QLineEdit *heightLineEdit;
  QComboBox *resolutionComboBox;
  QComboBox *mapPaletteComboBox;

  QAction *actionAntiAliasing;
  QAction *actionAitoff;
  QAction *actionFullSky;
  QAction *actionFillPixels;
  QAction *actionShowCoordinates;
  QAction *actionShowGrid;

  QGroupBox *pointsIOGroup;
  QPushButton *inputPointsButton;
  QLabel *currentPointsLabel;
  QComboBox *pointsCoordinateComboBox;
  QCheckBox *weightedPointsCheckBox;

  QGroupBox *pointsAppearanceGroup;
  QCheckBox *showPointsCheckBox;
  QCheckBox *filterPointsCheckBox;
  QGroupBox *pointsWeightGroup;
  QLineEdit *pointsWeightMinLineEdit;
  QLineEdit *pointsWeightMaxLineEdit;
  QGroupBox *pointsPaletteGroup;
  QComboBox *pointsPaletteComboBox;

  QLabel *dashLabel;
  QLabel *timesLabel;
  QErrorMessage *errorMessageDialog;
};

#endif
