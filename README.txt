CoroEval
========

CoroEval is a software for the evaluation of 3D coronary vessel reconstructions from clinical data. It runs on multiple operating systems and is designed to be independent of the imaging modality used. At this point, its purpose is the comparison of reconstruction algorithms or acquisition protocols, not the clinical diagnosis. Implemented metrics are vessel sharpness and diameter. All measurements are taken from the raw intensity data to be independent of display windowing functions. The user can either import a vessel centreline segmentation from other software, or perform a manual segmentation in CoroEval. An automated segmentation correction algorithm is provided to improve non-perfect centrelines.

The software is available in binary and source code form at http://www5.cs.fau.de/CoroEval/.

Authors
=======

Chris Schwemmer (chris.schwemmer@cs.fau.de)
Christoph Forman (christoph.forman@cs.fau.de)
Jens Wetzl (jens.wetzl@fau.de)

Citation
========

If you use CoroEval in your research, please cite its companion publication as follows:

CoroEval: a multi-platform, multi-modality tool for the evaluation of 3D coronary vessel reconstructions.
C Schwemmer, C Forman, J Wetzl, A Maier and J Hornegger.
Phys. Med. Biol. 59 (2014) 5163-5174.

Dependencies
============

CoroEval depends on these open source libraries:

* Qt (http://qt-project.org), tested with version 4.8.x
* Qwt (http://qwt.sourceforge.net), tested with version 6.1.0
* Grassroots DICOM library (GDCM, http://gdcm.sourceforge.net/wiki/index.php/Main_Page), tested with version 2.4.1
* libzip (http://www.nih.at/libzip/), tested with version 0.11.2
* Wild Magic 3 (http://geometrictools.com/). This version is no longer available online and therefore included in the source code package of CoroEval

Optionally, for mesh export:

* OpenMesh (http://openmesh.org), tested with version 3.2

Building from source
====================

CoroEval uses CMake (http://www.cmake.org) to generate build files for a wide variety of compiler environments, e.g. Visual Studio, Makefiles or Xcode. If you are unfamiliar with CMake, a brief introduction can be found at http://www.cmake.org/cmake/help/runningcmake.html. It is slightly outdated with respect to the supported types of compiler environments, but the workflow for creating build files is largely unchanged.

Once you have compiled the dependencies listed above or installed precompiled binaries, either include the root directories of all dependencies in the CMAKE_PREFIX_PATH environment variable or specify the *_DIR CMake variables when you run CMake to make sure the dependencies can be found by CoroEval. See the Find*.cmake scripts in the root directory of this software package for more details on how external dependencies are located.

By enabling the CMake option ENABLE_MESH_EXPORT, you will be able to use CoroEval to export meshes of segmented vessels. This requires OpenMesh in addition to the mandatory dependencies listed above.
