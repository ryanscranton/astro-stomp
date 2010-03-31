#include "stomp_viewer.h"
#include <QApplication>
#include <QMap>
#include <qdebug.h>

void usage() {
  qWarning() << "Usage: stomp_viewer [-SizeHint<color> <width>x<height>] ...";
  exit(1);
}

QMap<QString, QSize> parseCustomSizeHints(int argc, char **argv) {
  QMap<QString, QSize> result;

  for (int i=1;i<argc; ++i) {
    QString arg = QString::fromLocal8Bit(argv[i]);

    if (arg.startsWith(QLatin1String("-SizeHint"))) {
      QString name = arg.mid(9);
      if (name.isEmpty())
	usage();
      if (++i == argc)
	usage();
      QString sizeStr = QString::fromLocal8Bit(argv[i]);
      int idx = sizeStr.indexOf(QLatin1Char('x'));
      if (idx == -1)
	usage();
      bool ok;
      int w = sizeStr.left(idx).toInt(&ok);
      if (!ok)
	usage();
      int h = sizeStr.mid(idx + 1).toInt(&ok);
      if (!ok)
	usage();
      result[name] = QSize(w, h);
    }
  }

  return result;
 }

int main(int argc, char *argv[]) {
  QApplication app(argc, argv);
  QMap<QString, QSize> customSizeHints = parseCustomSizeHints(argc, argv);
  StompViewer stomp_viewer(customSizeHints,0);
  stomp_viewer.resize(1024, 680);
  stomp_viewer.show();
  return app.exec();
}

