QT += core
TEMPLATE = app

TARGET  = CoroEval

CONFIG += qt
CONFIG += qwt

INCLUDEPATH += src/
INCLUDEPATH += src/wm3

macx {
    QMAKE_LFLAGS += -F/System/Library/Frameworks/
    LIBS         += -framework CoreServices
    LIBS         += -framework Carbon
}
QMAKEFEATURES =  /usr/local/qwt-6.0.1/features

# Widget
HEADERS += src/imagewidget.h \
    src/Data.h \
    src/imagelabel.h \
    src/orthoviewer.h \
    src/wm3/Wm3Vector3.h \
    src/wm3/Wm3TTuple.h \
    src/wm3/Wm3TStack.h \
    src/wm3/Wm3TSet.h \
    src/wm3/Wm3TList.h \
    src/wm3/Wm3THashTable.h \
    src/wm3/Wm3THashSet.h \
    src/wm3/Wm3TArray.h \
    src/wm3/Wm3System.h \
    src/wm3/Wm3SingleCurve3.h \
    src/wm3/Wm3Segment3.h \
    src/wm3/Wm3Quaternion.h \
    src/wm3/Wm3Platforms.h \
    src/wm3/Wm3Memory.h \
    src/wm3/Wm3Matrix3.h \
    src/wm3/Wm3Math.h \
    src/wm3/Wm3Line3.h \
    src/wm3/Wm3Integrate1.h \
    src/wm3/Wm3FoundationLIB.h \
    src/wm3/Wm3Curve3.h \
    src/wm3/Wm3BSplineFitBasis.h \
    src/wm3/Wm3BSplineFit.h \
    src/wm3/Wm3BSplineCurve3.h \
    src/wm3/Wm3BSplineBasis.h \
    src/wm3/Wm3BandedMatrix.h \
    src/bsplineinterpolation.h \
    src/CImg.h \
    src/coronaryprofile.h \
    src/mainwindow.h \
    src/vesselsharpness.h \
    src/loaddialog.h \
    src/PlotWidget.h \
    src/Interval.h \
    src/sharpnessplot.h \
    src/settings.h \
    src/normalplot.h \
    src/pointwidget.h \
    src/measurementpoint.h \
    src/evaluationdialog.h \
    src/mipdialog.h
SOURCES += src/imagewidget.cpp \
    src/Data.cpp  \
    src/main.cpp \
    src/imagelabel.cpp \
    src/orthoviewer.cpp \
    src/wm3/Wm3Vector3.cpp \
    src/wm3/Wm3System.cpp \
    src/wm3/Wm3SingleCurve3.cpp \
    src/wm3/Wm3Quaternion.cpp \
    src/wm3/Wm3Memory.cpp \
    src/wm3/Wm3Matrix3.cpp \
    src/wm3/Wm3Math.cpp \
    src/wm3/Wm3Integrate1.cpp \
    src/wm3/Wm3Curve3.cpp \
    src/wm3/Wm3BSplineFitBasis.cpp \
    src/wm3/Wm3BSplineFit.cpp \
    src/wm3/Wm3BSplineCurve3.cpp \
    src/wm3/Wm3BSplineBasis.cpp \
    src/bsplineinterpolation.cpp \
    src/coronaryprofile.cpp \
    src/mainwindow.cpp \
    src/vesselsharpness.cpp \
    src/loaddialog.cpp \
    src/PlotWidget.cpp \
    src/sharpnessplot.cpp \
    src/settings.cpp \
    src/normalplot.cpp \
    src/pointwidget.cpp \
    src/measurementpoint.cpp \
    src/evaluationdialog.cpp \
    src/Interval.cpp \
    src/mipdialog.cpp
FORMS   +=      src/imagewidget.ui \
    src/orthoWidget.ui \
    src/coronaryprofile.ui \
    src/mainwindow.ui \
    src/loaddialog.ui \
    src/pointwidget.ui \
    src/evaluationdialog.ui

RESOURCES += \
    src/resource.qrc
