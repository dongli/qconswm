cmake_minimum_required (VERSION 2.8)

project (qconswm Fortran)

find_package (MPI)
if (NOT ${MPI_Fortran_FOUND})
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DNO_MPI")
else ()
  set (CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif ()

set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")

find_package(PkgConfig)
if (${PKG_CONFIG_FOUND})
  pkg_search_module (DATETIME REQUIRED datetime-fortran)
  if (${DATETIME_FOUND})
    include_directories (${DATETIME_INCLUDE_DIRS})
    link_directories (${DATETIME_LIBRARY_DIRS})
  else ()
    message (FATAL_ERROR "Pkg-config could not find datetime-fortran library!")
  endif ()
else ()
  message (FATAL_ERROR "Unable to find datetime-fortran library!")
endif ()

set (source_files
  "src/params_mod.F90"
  "src/mesh_mod.F90"
  "src/parallel_mod.F90"
  "src/dycore_mod.F90"
  "src/qconswm.F90"
)

add_executable (qconswm.exe ${source_files})
target_link_libraries (qconswm.exe ${DATETIME_LIBRARIES})
