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
find_package(Qt5Widgets)
if (Qt5Widgets_FOUND)
    set(QT_LIBRARIES Qt5::Widgets)
    set(USE_QT5 ON)
    find_package(Qwt REQUIRED)
    INCLUDE_DIRECTORIES(${QWT_INCLUDE_DIRS})
else()
    set(USE_QT5 OFF)
    find_package(Qt4 REQUIRED)
    include(${QT_USE_FILE})
    add_definitions(${QT_DEFINITIONS})
    FIND_PACKAGE(Qwt6 REQUIRED)
endif()


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
    ${CoroEval_SOURCE_DIR}/wm3
)

IF(WIN32)
    ADD_DEFINITIONS(-DQT_DLL)
    ADD_DEFINITIONS(-DQWT_DLL)
ENDIF(WIN32)


# Source setup


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


set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wno-shift-negative-value -std=c++14")

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
    ${CoroEval_RCS}
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
    ${QWT_LIBRARIES}
    wm3
)

target_link_libraries(${PROJECT_NAME} ${QWT_LIBRARIES} wm3)

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
        gdcmMSFF
        gdcmopenjp2
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

add_subdirectory(wm3)
