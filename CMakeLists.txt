cmake_minimum_required (VERSION 2.8.8)

project (korcgpu_proj Fortran)

#  Set directories to build all binary files.
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build/bin)

#  Fortran specific settings.
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/modules)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/lib)

#  Set definitions.
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
   add_definitions (-DDBG_CHECK)
endif()

#  Configure OpenMP
find_package (OpenMP)
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_Fortran_FLAGS}")

#  Configure PSPLINE
option (USE_PSPLINE "Use PSPLINE interpolation" OFF)
if (${USE_PSPLINE})
   include_directories (${PSPLINE_INCLUDE_PATH})
   set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PSPLINE_LINK_FLAGS}")
   add_definitions (-DPSPLINE)
endif ()

add_subdirectory (src)
