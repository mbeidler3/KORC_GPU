cmake_minimum_required (VERSION 2.8.8)

project (korcgpu_proj Fortran)

#  Set directories to build all binary files.
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build/bin)

#  Fortran specific settings.
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/modules)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/lib)

add_definitions (-DDOUBLE_PRECISION)

#  Set definitions.
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
   add_definitions (-DDBG_CHECK)
endif()

#Configure OpenMP
option (USE_OMP "Use OpenOMP" OFF)
if (${USE_OMP})
   find_package (OpenMP)
   set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
   set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_Fortran_FLAGS}")
   add_definitions (-DOMP)
endif()

option (USE_ACC "Use OpenACC" OFF)
if (${USE_ACC})
   add_definitions (-DACC)
endif()

#  Configure HDF5
#find_package (HDF5 REQUIRED COMPONENTS Fortran)
#include_directories (${HDF5_INCLUDE_DIRS})

#  Configure MPI
#find_package (MPI REQUIRED)
#set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
#include_directories (${MPI_Fortran_INCLUDE_PATH})

#  Configure PSPLINE
option (USE_PSPLINE "Use PSPLINE interpolation" OFF)
if (${USE_PSPLINE})
   #include_directories (${PSPLINE_INCLUDE_PATH})
   #set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PSPLINE_LINK_FLAGS}")
   add_definitions (-DPSPLINE)
endif ()

add_subdirectory (src)
