#ifndef CONFIG_H
#define CONFIG_H

#define CMAKE_Fortran_COMPILER_ID GNU
#define CMAKE_Fortran_FLAGS -fconvert=big-endian -fno-realloc-lhs -fno-second-underscore  -pipe -fallow-argument-mismatch -fdefault-real-8 -fdefault-double-8 
#define CMAKE_Fortran_BUILD_TYPE_FLAGS -O3
#define Fortran_COMPILER_NAME gfortran
#define Fortran_VERSION 12.1.0
#define HAVE_BLAS
#define HAVE_CMAKE
#define HAVE_FC_HDF5
#define HAVE_HDF5
#define HAVE_LAPACK
#define HAVE_NIMBND
#define HAVE_SUPERLU_SEQ
#define HAVE_OCULUS
#define HAVE_BSC
#define HOSTTYPE Darwin-arm64
#define PROJECT_REV unknown
#define PROJECT_URL unknown
#define USE_LE_SURFACE
#define UQHOSTNAME MAC130285
#define __gfortran

#endif 
