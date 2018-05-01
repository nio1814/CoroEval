find_path(OpenMesh_INCLUDE_DIR OpenMesh HINTS /usr/include /usr/local/include /opt/local/include PATH_SUFFIXES OpenMesh)
find_library(OpenMesh_LIBRARY1 NAMES OpenMeshCore HINTS /usr/lib /usr/local/lib /opt/local/lib PATH_SUFFIXES OpenMesh)
find_library(OpenMesh_LIBRARY2 NAMES OpenMeshTools HINTS /usr/lib /usr/local/lib /opt/local/lib PATH_SUFFIXES OpenMesh)

set(OpenMesh_INCLUDE_DIRS ${OpenMesh_INCLUDE_DIR})
set(OpenMesh_LIBRARIES ${OpenMesh_LIBRARY1} ${OpenMesh_LIBRARY2})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenMesh DEFAULT_MSG OpenMesh_INCLUDE_DIR OpenMesh_LIBRARY1 OpenMesh_LIBRARY2)

mark_as_advanced(OpenMesh_INCLUDE_DIR OpenMesh_LIBRARY1 OpenMesh_LIBRARY2)
