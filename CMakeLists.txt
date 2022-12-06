cmake_minimum_required (VERSION 2.8.8)

project (korc_proj Fortran)

#  Set directories to build all binary files.
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build/bin)

#  Fortran specific settings.
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/modules)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/lib)

#  Set definitions.
add_definitions (-DDOUBLE_PRECISION)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
   add_definitions (-DDBG_CHECK)
endif()

#  Configure OpenMP
find_package (OpenMP)
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_Fortran_FLAGS}")

#  Configure MPI
find_package (MPI REQUIRED)
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
include_directories (${MPI_Fortran_INCLUDE_PATH})

#  Configure HDF5
set( HDF5_USE_STATIC_LIBRARIES ON )
find_package( HDF5 COMPONENTS Fortran REQUIRED )
include_directories( SYSTEM ${HDF5_INCLUDE_DIRS} )

add_subdirectory (src)
