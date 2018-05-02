# Find the Qwt 6.x includes and library, linked to Qt4
#
# On Windows it makes these assumptions:
#    - the Qwt DLL is where the other DLLs for Qt are (QT_DIR\bin) or in the path
#    - the Qwt .h files are in QT_DIR\include\Qwt or in the path
#    - the Qwt .lib is where the other LIBs for Qt are (QT_DIR\lib) or in the path
#
# Qwt_INCLUDE_DIR - where to find qwt.h if Qwt
# Qwt-Qt4_LIBRARY - The Qwt6 library linked against Qt4 (if it exists)
# Qwt-Qt4_FOUND   - Qwt6 was found and uses Qt4
# Qwt_FOUND - Set to TRUE if Qwt6 was found (linked to Qt4)

# Copyright (c) 2007, Pau Garcia i Quiles, <pgquiles@elpauer.org>
# Changes by Chris Schwemmer, Universität Erlangen-Nürnberg, <chris.schwemmer@cs.fau.de>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

# Condition is "(A OR B) AND C", CMake does not support parentheses but it evaluates left to right
if (NOT USE_QT5)
    IF(Qwt-Qt4_LIBRARY AND Qwt_INCLUDE_DIR)
        SET(Qwt_FIND_QUIETLY TRUE)
    ENDIF(Qwt-Qt4_LIBRARY AND Qwt_INCLUDE_DIR)

    FIND_PACKAGE( Qt4 REQUIRED QUIET )

    IF( QT4_FOUND )
      # Is Qwt installed? Look for header files
      FIND_PATH( Qwt_INCLUDE_DIR qwt.h PATHS ${QT_INCLUDE_DIR} /usr/local/Cellar/qwt/6.1.0/lib/qwt.framework/Versions/6/Headers PATH_SUFFIXES qwt qwt6 qwt-qt4 qwt6-qt4 include qwt/include qwt6/include qwt-qt4/include qwt6-qt4/include ENV PATH)
endif(NOT USE_QT5)
      
	# Find Qwt version
	IF( Qwt_INCLUDE_DIR )
		FILE( READ ${Qwt_INCLUDE_DIR}/qwt_global.h QWT_GLOBAL_H )
		STRING( REGEX MATCH "#define *QWT_VERSION *(0x06*)" QWT_IS_VERSION_6 ${QWT_GLOBAL_H})
		
		IF( QWT_IS_VERSION_6 )
		STRING(REGEX REPLACE ".*#define[\\t\\ ]+QWT_VERSION_STR[\\t\\ ]+\"([0-9]+\\.[0-9]+\\.[0-9]+)\".*" "\\1" Qwt_VERSION "${QWT_GLOBAL_H}")

		# Find Qwt library linked to Qt4
		FIND_LIBRARY( Qwt-Qt4_TENTATIVE_LIBRARY NAMES qwt6-qt4 qwt-qt4 qwt6 qwt HINTS ${Qwt_INCLUDE_DIR}/../lib )
		FIND_LIBRARY( Qwt-Qt4_TENTATIVE_LIBRARY_DEBUG NAMES qwt6-qt4d qwt-qt4d qwt6d qwtd HINTS ${Qwt_INCLUDE_DIR}/../lib )
		IF( UNIX AND NOT CYGWIN AND NOT APPLE)
			IF( Qwt-Qt4_TENTATIVE_LIBRARY )
				#MESSAGE("Qwt-Qt4_TENTATIVE_LIBRARY = ${Qwt-Qt4_TENTATIVE_LIBRARY}")
				EXECUTE_PROCESS( COMMAND "ldd" ${Qwt-Qt4_TENTATIVE_LIBRARY} OUTPUT_VARIABLE Qwt-Qt4_LIBRARIES_LINKED_TO )
				STRING( REGEX MATCH "QtCore" Qwt_IS_LINKED_TO_Qt4 ${Qwt-Qt4_LIBRARIES_LINKED_TO})
				IF( Qwt_IS_LINKED_TO_Qt4 )
					SET( Qwt-Qt4_LIBRARY ${Qwt-Qt4_TENTATIVE_LIBRARY} )
					SET( Qwt-Qt4_FOUND TRUE )
					IF (NOT Qwt_FIND_QUIETLY)
						MESSAGE( STATUS "Found Qwt version ${Qwt_VERSION} linked to Qt4" )
					ENDIF (NOT Qwt_FIND_QUIETLY)
				ENDIF( Qwt_IS_LINKED_TO_Qt4 )
			ENDIF( Qwt-Qt4_TENTATIVE_LIBRARY )
		ELSEIF(APPLE)
			SET(Qwt-Qt4_LIBRARY ${Qwt-Qt4_TENTATIVE_LIBRARY})
			SET(Qwt-Qt4_FOUND TRUE )
			IF (NOT Qwt_FIND_QUIETLY)
				MESSAGE( STATUS "Found Qwt version ${Qwt_VERSION} linked to Qt4" )
			ENDIF (NOT Qwt_FIND_QUIETLY)
		ELSE(APPLE)
		# Assumes qwt.dll is in the Qt dir
			SET( Qwt-Qt4_LIBRARY optimized ${Qwt-Qt4_TENTATIVE_LIBRARY} debug ${Qwt-Qt4_TENTATIVE_LIBRARY_DEBUG})
			SET( Qwt-Qt4_FOUND TRUE )
			IF (NOT Qwt_FIND_QUIETLY)
				MESSAGE( STATUS "Found Qwt version ${Qwt_VERSION} linked to Qt4" )
			ENDIF (NOT Qwt_FIND_QUIETLY)
		ENDIF( UNIX AND NOT CYGWIN AND NOT APPLE)
			
		ENDIF( QWT_IS_VERSION_6 )
		
		IF( Qwt-Qt4_FOUND )
			SET( Qwt_FOUND TRUE )
		ENDIF( Qwt-Qt4_FOUND )
		
		MARK_AS_ADVANCED( Qwt_INCLUDE_DIR Qwt-Qt4_LIBRARY )
	ENDIF( Qwt_INCLUDE_DIR )

   	IF (NOT Qwt_FOUND AND Qwt_FIND_REQUIRED)
      		MESSAGE(FATAL_ERROR "Could not find Qwt 6.x")
   	ENDIF (NOT Qwt_FOUND AND Qwt_FIND_REQUIRED)

ENDIF( QT4_FOUND )
