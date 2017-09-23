cmake_minimum_required (VERSION 2.8)

project (qconswm Fortran)

set (CMAKE_Fortran_FLAGS "-fdefault-real-8")

set (source_files
  "src/params_mod.F90"
  "src/mesh_mod.F90"
  "src/dycore_mod.F90"
  "src/qconswm.F90"
)

add_executable (qconswm.exe ${source_files})
