# CoroEval
#
# An evaluation tool for coronary artery reconstructions.
#
# Copyright © 2014:
#
# Christoph Forman, Universität Erlangen-Nürnberg
# christoph.forman@cs.fau.de
#
# Chris Schwemmer, Universität Erlangen-Nürnberg
# chris.schwemmer@cs.fau.de
#
# Jens Wetzl, Universität Erlangen-Nürnberg
# jens.wetzl@cs.fau.de
#
# CoroEval is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CoroEval is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CoroEval.  If not, see <http://www.gnu.org/licenses/>.

PROJECT(CoroEval)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CMAKE_MODULE_PATH ${CoroEval_SOURCE_DIR})

# Requirements

FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

FIND_PACKAGE(Qwt6 REQUIRED)

FIND_PACKAGE(LibZip REQUIRED)

FIND_PACKAGE(GDCM REQUIRED)
INCLUDE(${GDCM_USE_FILE})

OPTION(ENABLE_MESH_EXPORT "Enable exporting meshes from vessel segmentations (requires OpenMesh)" OFF)

IF(ENABLE_MESH_EXPORT)
    FIND_PACKAGE(OpenMesh REQUIRED)

    ADD_DEFINITIONS(-DENABLE_MESH_EXPORT)
    INCLUDE_DIRECTORIES(${OpenMesh_INCLUDE_DIRS})
    ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
ENDIF(ENABLE_MESH_EXPORT)

#TODO: Frameworks for MacOS X

# Compiler and linker environment

INCLUDE_DIRECTORIES(
    ${QT_INCLUDE_DIR}
    ${QT_INCLUDE_DIR}/Qt
    ${Qwt_INCLUDE_DIR}
    ${ZLIB_INCLUDE_DIR}
    ${libzip_INCLUDE_DIRS}
    ${CoroEval_BINARY_DIR}
    ${CoroEval_SOURCE_DIR}/src
    ${CoroEval_SOURCE_DIR}/src/wm3
)

IF(WIN32)
    ADD_DEFINITIONS(-DQT_DLL)
    ADD_DEFINITIONS(-DQWT_DLL)
ENDIF(WIN32)


# Source setup

SET(WM3_HDRS
    src/wm3/Wm3Vector2.h
    src/wm3/Wm3Vector3.h
    src/wm3/Wm3Vector4.h
    src/wm3/Wm3GVector.h
    src/wm3/Wm3TTuple.h
    src/wm3/Wm3TStack.h
    src/wm3/Wm3TSet.h
    src/wm3/Wm3TList.h
    src/wm3/Wm3THashTable.h
    src/wm3/Wm3THashSet.h
    src/wm3/Wm3TArray.h
    src/wm3/Wm3System.h
    src/wm3/Wm3SingleCurve3.h
    src/wm3/Wm3Segment3.h
    src/wm3/Wm3Quaternion.h
    src/wm3/Wm3Platforms.h
    src/wm3/Wm3Memory.h
    src/wm3/Wm3Matrix2.h
    src/wm3/Wm3Matrix3.h
    src/wm3/Wm3Matrix4.h
    src/wm3/Wm3GMatrix.h
    src/wm3/Wm3Math.h
    src/wm3/Wm3Line3.h
    src/wm3/Wm3Integrate1.h
    src/wm3/Wm3FoundationLIB.h
    src/wm3/Wm3Curve3.h
    src/wm3/Wm3BSplineFitBasis.h
    src/wm3/Wm3BSplineFit.h
    src/wm3/Wm3BSplineCurve3.h
    src/wm3/Wm3BSplineBasis.h
    src/wm3/Wm3BandedMatrix.h
    src/wm3/Wm3Eigen.h
	src/wm3/Wm3NoniterativeEigen3x3.h
)

SET(WM3_SRCS
    src/wm3/Wm3Vector2.cpp
    src/wm3/Wm3Vector3.cpp
    src/wm3/Wm3Vector4.cpp
    src/wm3/Wm3System.cpp
    src/wm3/Wm3SingleCurve3.cpp
    src/wm3/Wm3Quaternion.cpp
    src/wm3/Wm3Memory.cpp
    src/wm3/Wm3Matrix2.cpp
    src/wm3/Wm3Matrix3.cpp
    src/wm3/Wm3Matrix4.cpp
    src/wm3/Wm3Math.cpp
    src/wm3/Wm3Integrate1.cpp
    src/wm3/Wm3Curve3.cpp
    src/wm3/Wm3BSplineFitBasis.cpp
    src/wm3/Wm3BSplineFit.cpp
    src/wm3/Wm3BSplineCurve3.cpp
    src/wm3/Wm3BSplineBasis.cpp
    src/wm3/Wm3Eigen.cpp
	src/wm3/Wm3NoniterativeEigen3x3.cpp
)

SET(CoroEval_HDRS
    src/imagewidget.h
    src/Data.h
    src/imagelabel.h
    src/orthoviewer.h
    src/bsplineinterpolation.h
    src/CImg.h
    src/coronaryprofile.h
	src/coronaryprofile2.h
    src/mainwindow.h
    src/vesselsharpness.h
    src/loaddialog.h
    src/PlotWidget.h
    src/Interval.h
    src/sharpnessplot.h
    src/settings.h
    src/normalplot.h
    src/pointwidget.h
    src/measurementpoint.h
    src/evaluationdialog.h
    src/mipdialog.h
    src/truncatedialog.h
    src/settingsdialog.h
	src/ZipFile.h
	src/MPFilter.h
)

SET(CoroEval_SRCS
    src/imagewidget.cpp
    src/Data.cpp
    src/main.cpp
    src/imagelabel.cpp
    src/orthoviewer.cpp
    src/bsplineinterpolation.cpp
    src/coronaryprofile.cpp
	src/coronaryprofile2.cpp
    src/mainwindow.cpp
    src/vesselsharpness.cpp
    src/loaddialog.cpp
    src/PlotWidget.cpp
    src/sharpnessplot.cpp
    src/settings.cpp
    src/normalplot.cpp
    src/pointwidget.cpp
    src/measurementpoint.cpp
    src/evaluationdialog.cpp
    src/Interval.cpp
    src/mipdialog.cpp
    src/truncatedialog.cpp
    src/settingsdialog.cpp
	src/ZipFile.cpp
	src/MPFilter.cpp
)

IF(ENABLE_MESH_EXPORT)
	SET(CoroEval_HDRS
		${CoroEval_HDRS}
		src/Triangulation.h
	)
	SET(CoroEval_SRCS
		${CoroEval_SRCS}
		src/Triangulation.cpp
	)
ENDIF(ENABLE_MESH_EXPORT)

QT4_WRAP_UI(CoroEval_UIS_H
    src/imagewidget.ui
    src/orthoWidget.ui
    src/coronaryprofile.ui
	src/coronaryprofile2.ui
    src/mainwindow.ui
    src/loaddialog.ui
    src/pointwidget.ui
    src/evaluationdialog.ui
    src/truncatedialog.ui
    src/settingsdialog.ui
)

SET(CoroEval_MOC_SRCS
    src/coronaryprofile.h
	src/coronaryprofile2.h
	src/imagelabel.h
    src/imagewidget.h
    src/loaddialog.h
    src/mainwindow.h
    src/normalplot.h
    src/orthoviewer.h
    src/PlotWidget.h
    src/sharpnessplot.h
    src/pointwidget.h
    src/evaluationdialog.h
    src/mipdialog.h
    src/truncatedialog.h
    src/settingsdialog.h
)


SET(CoroEval_RCS
    src/resource.qrc
)

SET(CoroPlot_HDRS
    src/Data.h
    src/CImg.h
    src/bsplineinterpolation.h
    src/vesselsharpness.h
    src/measurementpoint.h
    src/Interval.h
    src/settings.h
    src/normalplot.h
    src/PlotWidget.h
    src/evaldialog.h
)

SET(CoroPlot_SRCS
    src/Data.cpp
    src/imagelabel.cpp
    src/mainCoroPlot.cpp
    src/bsplineinterpolation.cpp
    src/vesselsharpness.cpp
    src/measurementpoint.cpp
    src/Interval.cpp
    src/settings.cpp
    src/normalplot.cpp
    src/PlotWidget.cpp
    src/evaldialog.cpp
)

SET(CoroPlot_MOC_SRCS
    src/imagelabel.h
    src/normalplot.h
    src/PlotWidget.h
    src/evaldialog.h
)

QT4_WRAP_CPP(MOC_CXX
    ${CoroEval_MOC_SRCS}
)

QT4_WRAP_CPP(MOC_CXX2
    ${CoroPlot_MOC_SRCS}
)

QT4_ADD_RESOURCES(CoroEval_RCS_CXX
    ${CoroEval_RCS}
)

IF(WIN32)
ADD_EXECUTABLE(CoroEval
    ${WM3_SRCS}
    ${WM3_HDRS}
    ${CoroEval_SRCS}
    ${CoroEval_HDRS}
    ${CoroEval_UIS_H}
    ${MOC_CXX}
    ${CoroEval_RCS_CXX}
	${CoroEval_SOURCE_DIR}/src/icon.rc
)
ELSE(WIN32)
ADD_EXECUTABLE(CoroEval
    ${WM3_SRCS}
    ${WM3_HDRS}
    ${CoroEval_SRCS}
    ${CoroEval_HDRS}
    ${CoroEval_UIS_H}
    ${MOC_CXX}
    ${CoroEval_RCS_CXX}
)
ENDIF(WIN32)

ADD_EXECUTABLE(CoroPlot
    ${WM3_SRCS}
    ${WM3_HDRS}
    ${CoroPlot_SRCS}
    ${CoroPlot_HDRS}
    ${MOC_CXX2}
    ${CoroEval_RCS_CXX}
)

TARGET_LINK_LIBRARIES(CoroPlot
    ${QT_LIBRARIES}
    ${Qwt-Qt4_LIBRARY}
)

IF(WIN32)
    TARGET_LINK_LIBRARIES(CoroEval
        gdcmCommon
        gdcmDICT
        gdcmDSED
        gdcmIOD
        gdcmMSFF
        gdcmzlib
        gdcmexpat
        gdcmjpeg8
        gdcmjpeg12
        gdcmjpeg16
        gdcmopenjpeg
        gdcmcharls
        ${QT_LIBRARIES}
        ${Qwt-Qt4_LIBRARY}
        ${libzip_LIBRARY}
    )
ELSE(WIN32)
    TARGET_LINK_LIBRARIES(CoroEval
        gdcmCommon
        gdcmDICT
        gdcmDSED
        gdcmIOD
        gdcmMSFF
        gdcmzlib
        gdcmexpat
        gdcmuuid
        gdcmjpeg8
        gdcmjpeg12
        gdcmjpeg16
        gdcmopenjpeg
        gdcmcharls
        ${QT_LIBRARIES}
        ${Qwt-Qt4_LIBRARY}
        ${libzip_LIBRARY}
    )
ENDIF(WIN32)

IF(ENABLE_MESH_EXPORT)
	TARGET_LINK_LIBRARIES(CoroEval
		${OpenMesh_LIBRARIES}
	)
ENDIF(ENABLE_MESH_EXPORT)

IF(APPLE)
	FIND_LIBRARY(CORE_SERVICES_LIBRARY CoreServices)
	FIND_LIBRARY(CARBON_LIBRARY Carbon)
	TARGET_LINK_LIBRARIES(CoroEval ${CORE_SERVICES_LIBRARY} ${CARBON_LIBRARY})
	TARGET_LINK_LIBRARIES(CoroPlot ${CORE_SERVICES_LIBRARY} ${CARBON_LIBRARY})
ENDIF(APPLE)
