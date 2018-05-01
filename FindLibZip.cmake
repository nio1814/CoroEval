# Find the libzip includes and library
#
# libzip_INCLUDE_DIRS - where to find zip.h
# libzip_LIBRARY - The libzip library
# libzip_FOUND   - libzip was found

IF(libzip_LIBRARY AND libzip_INCLUDE_DIRS)
	SET(libzip_FIND_QUIETLY TRUE)
ENDIF(libzip_LIBRARY AND libzip_INCLUDE_DIRS)

FIND_PACKAGE(ZLIB REQUIRED QUIET)

# Is libzip installed? Look for header files
FIND_PATH(libzip_INCLUDE_DIR zip.h PATHS PATH_SUFFIXES libzip zip include libzip/include zip/include ENV PATH)

IF(libzip_INCLUDE_DIR )
	FIND_LIBRARY(libzip_LIBRARY NAMES zip HINTS ${libzip_INCLUDE_DIR}/../lib ${libzip_INCLUDE_DIR}/../${CMAKE_BUILD_TYPE})
	SET(libzip_FOUND TRUE )
	IF (NOT libzip_FIND_QUIETLY)
		MESSAGE(STATUS "Found libzip")
	ENDIF (NOT libzip_FIND_QUIETLY)

	IF(MSVC)
		ADD_DEFINITIONS("-D_CRT_SECURE_NO_WARNINGS")
		ADD_DEFINITIONS("-DZLIB_WINAPI")
	ENDIF(MSVC)

	GET_FILENAME_COMPONENT(_libzip_LIBPATH ${libzip_LIBRARY} PATH)
	SET(libzip_INCLUDE_DIRS "${libzip_INCLUDE_DIR};${_libzip_LIBPATH}/libzip/include" CACHE PATH "")

	MARK_AS_ADVANCED(libzip_INCLUDE_DIR libzip_INCLUDE_DIRS libzip_LIBRARY)
ENDIF(libzip_INCLUDE_DIR)

IF (NOT libzip_FOUND AND libzip_FIND_REQUIRED)
	MESSAGE(FATAL_ERROR "Could not find libzip")
ENDIF (NOT libzip_FOUND AND libzip_FIND_REQUIRED)

