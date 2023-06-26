#!/bin/bash

_error() {

  local readonly msg="$1"

  echo "."
  echo "."
  echo "launchkit >> [Error] << ${msg}"
  echo "."
  echo "."

  return 1
}


_error_fatal() {

  local readonly msg="$1"

  echo "."
  echo "."
  echo "launchkit >> [Fatal Error] << ${msg}"
  echo "."
  echo "."

  exit 1
}

_platform_korcgpu_match() {

  [[ $HOSTNAME =~ korcgpu ]] || return 1
}

_platform_korcgpu_setup() {
  export CC=nvcc
  export CXX=nvc++
  export FC=nvfortran
}

_platform_mattlaptop_match() {

  [[ $(uname -n) =~ MBP115573 ]] || return 1
}


_platform_mattlaptop_setup() {
  export CC=/usr/local/bin/gcc-13
  export CXX=/usr/local/bin/g++-13
  export FC=/usr/local/bin/gfortran-13
}


_platform_perlmutter_match() {

  [[ $LMOD_SYSTEM_NAME =~ perlmutter ]] || return 1
}


_platform_perlmutter_setup() {

  module load cmake
  module load nvidia
}

_scenario_gpu_trial1_korcgpu_build() {

  BUILD_TYPE=Release
  HARDWARE_TYPE=GPU
 
  echo cmake configured to generate a $BUILD_TYPE build.

  rm -f CMakeCache.txt
  rm -rf CMakeFiles

  rm ./src/*.mod

  rm -rf ./build
  mkdir ./build && cd $_

  if [ "$HARDWARE_TYPE" == CPU ]
  then
    cmake \
        -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
        -DUSE_OMP=OFF \
        -DCMAKE_Fortran_FLAGS="-Mfree -Mpreprocess -r8 -fPIC -O3 -mp" \
        -DCMAKE_Fortran_FLAGS_DEBUG="-g -Minfo=all -Minstrument -traceback -lnvhpcwrapnvtx" \
        -DUSE_PSPLINE=ON \
      ..
  elif [ "$HARDWARE_TYPE" == GPU ]
  then
    cmake \
        -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
        -DUSE_OMP=OFF \
        -DUSE_ACC=ON \
        -DCMAKE_Fortran_FLAGS="-fast -acc=gpu -gpu=deepcopy,cc70 -Mfree -fPIC -O3 -Mpreprocess" \
        -DCMAKE_Fortran_FLAGS_DEBUG="-g -Minfo=all -Minline -Minstrument -traceback -lnvhpcwrapnvtx" \
        -DUSE_PSPLINE=ON \
      ..
  else
    echo " choose a valid hardware type "
  fi

  [[ $? -eq 0 ]] && make VERBOSE=1

  LAUNCHKIT_RUNTIME_PATH_EXE="$(pwd -P)/build/bin/xkorcgpu"
  cd - >/dev/null
}

_scenario_gpu_trial1_mattlaptop_build() {

  BUILD_TYPE=Debug
 
  echo cmake configured to generate a $BUILD_TYPE build.

  rm -f CMakeCache.txt
  rm -rf CMakeFiles

  rm ./src/*.mod

  rm -rf ./build
  mkdir ./build && cd $_

  cmake \
      -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
      -DUSE_OMP=OFF \
      -DCMAKE_Fortran_FLAGS="-O3 -fopenmp -fallow-argument-mismatch -malign-double" \
      -DCMAKE_C_FLAGS="-O3 -fopenmp -malign-double"  \
      -DCMAKE_CXX_FLAGS="-O3 -fopenmp -malign-double" \
      -DCMAKE_Fortran_FLAGS_DEBUG="-g -ffpe-trap=invalid,zero,overflow -fbacktrace" \
      -DCMAKE_C_FLAGS_DEBUG="-g -g3" \
      -DCMAKE_CXX_FLAGS_DEBUG="-g -g3" \
      -DUSE_PSPLINE=ON \
    ..

  [[ $? -eq 0 ]] && make VERBOSE=1

  LAUNCHKIT_RUNTIME_PATH_EXE="$(pwd -P)/build/bin/xkorcgpu"
  cd - >/dev/null
}


_scenario_gpu_trial1_mattlaptop_run() {

  echo
}


_scenario_gpu_trial1_perlmutter_build() {

  BUILD_TYPE=Debug
  HARDWARE_TYPE=GPU
 
  echo cmake configured to generate a $BUILD_TYPE build.

  rm -f CMakeCache.txt
  rm -rf CMakeFiles

  rm ./src/*.mod

  rm -rf ./build
  mkdir ./build && cd $_

  if [ "$HARDWARE_TYPE" == CPU ]
  then
    cmake \
        -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
        -DUSE_OMP=ON \
        -DCMAKE_Fortran_FLAGS="-Mfree -Mpreprocess -r8 -fPIC -O3 -mp" \
        -DCMAKE_Fortran_FLAGS_DEBUG="-g -Minfo=all -Minstrument -traceback -lnvhpcwrapnvtx" \
        -DUSE_PSPLINE=ON \
        -DPSPLINE_INCLUDE_PATH=/global/cfs/cdirs/m3236/build_unstable/PSPLINE_GH/pspline/build/include \
        -DPSPLINE_LINK_FLAGS="-L/global/cfs/cdirs/m3236/build_unstable/PSPLINE_GH/pspline/build/lib -lpspline" \
        -DPSPLINE_LIBRARIES=/global/cfs/cdirs/m3236/build_unstable/PSPLINE_GH/pspline/build/lib/libpspline.a   \
      ..
  elif [ "$HARDWARE_TYPE" == GPU ]
  then
    cmake \
        -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
        -DUSE_OMP=OFF \
        -DUSE_ACC=ON \
        -DCMAKE_Fortran_FLAGS="-fast -acc=gpu -gpu=cc80 -Mfree -fPIC -O3 -Mpreprocess" \
        -DCMAKE_Fortran_FLAGS_DEBUG="-g -Minfo=all -Minstrument -traceback -lnvhpcwrapnvtx" \
        -DUSE_PSPLINE=ON \
        -DPSPLINE_INCLUDE_PATH=/global/cfs/cdirs/m3236/build_gpu_test/PSPLINE_GH/pspline/build/include \
        -DPSPLINE_LINK_FLAGS="-L/global/cfs/cdirs/m3236/build_gpu_test/PSPLINE_GH/pspline/build/lib -lpspline" \
        -DPSPLINE_LIBRARIES=/global/cfs/cdirs/m3236/build_gpu_test/PSPLINE_GH/pspline/build/lib/libpspline.a   \
      ..
  else
    echo " choose a valid hardware type "
  fi

  [[ $? -eq 0 ]] && make VERBOSE=1

  LAUNCHKIT_RUNTIME_PATH_EXE="$(pwd -P)/build/bin/xkorcgpu"
  cd - >/dev/null
}


_scenario_gpu_trial1_perlmutter_run() {

  cd $(mktemp -d) || return 1
  mkdir .launchkit && cd $_
  launchkit_dugout_dir="$(pwd -P)"

  cat > ./slurm_job.sh <<-__EOF__
	#!/bin/bash

	#SBATCH -A ntrain1
	#SBATCH -C gpu
	#SBATCH -q regular
	#SBATCH -t 1:00:00
	#SBATCH -N 1
	#SBATCH --ntasks-per-node=1
	#SBATCH -c 128
	#SBATCH --gpus-per-task=1
	#SBATCH --gpu-bind=none

	export SLURM_CPU_BIND="cores"
	srun nsys profile -o thinger --stats=true $LAUNCHKIT_RUNTIME_PATH_EXE

	__EOF__

  echo "."; echo "."; cat ./slurm_job.sh
  echo "."; echo "."

  sbatch ./slurm_job.sh

  cd - >/dev/null

  # cleanup temp dir
  [[ ${launchkit_dugout_dir%/*} =~ .+/tmp\..+ ]] && rm -rf ${launchkit_dugout_dir%/*}
}


__platform_match() {

  local match_count=0

  [[ -n $LAUNCHKIT_PLATFORM ]] && return

  _platform_perlmutter_match && let "match_count+=1" && LAUNCHKIT_PLATFORM="perlmutter"
  _platform_mattlaptop_match && let "match_count+=1" && LAUNCHKIT_PLATFORM="mattlaptop"
  _platform_korcgpu_match && let "match_count+=1" && LAUNCHKIT_PLATFORM="korcgpu"

  [[ ${match_count} -eq 1 ]] || _error_fatal "unable to select platform"
}


__platform_setup() {

  local readonly platform="$LAUNCHKIT_PLATFORM"

  _platform_${platform}_setup
}


__build_and_run() {

  local readonly platform="$LAUNCHKIT_PLATFORM"

  _scenario_gpu_trial1_${platform}_build &> compiler_report.txt
  #_scenario_gpu_trial1_${platform}_run
}


LAUNCHKIT_RUNTIME_PATH_EXE=""
LAUNCHKIT_PLATFORM=""


__platform_match
__platform_setup

__build_and_run



