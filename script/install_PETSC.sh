#! /bin/bash

# check if the directory $1/PETSC exist

CXX=$4
CC=$3
F77=$5
FC=$6

if [ -d "$1/PETSC" -a -f "$1/PETSC/include/petsc.h" ]; then
  echo "PETSC is already installed"
  exit 0
fi

# Detect gcc pr clang

source script/discover_os
source script/solve_python
discover_os

function haveProg() {
    [ -x "$(command -v $1)" ]
}

if haveProg python2; then
  python_command=python2
else
  python_command=python
fi


##### if we are on osx we use gsed

ldflags_petsc=
if [ x"$platform" == x"osx" ]; then
  sed_command=gsed
  ldflags_petsc=
else
  sed_command=sed
  ldflags_petsc=
fi


####

## If some dependencies has been installed feed them to PETSC

MUMPS_extra_libs=""

configure_options=""
configure_options_superlu=""
configure_trilinos_options=" -D TPL_ENABLE_MPI=ON "
configure_options_hypre=""

### Here we install OpenBLAS and SUITESPARSE

configure_options="$configure_options --download-metis --download-parmetis"

if [ -d "$1/BOOST" ]; then
  configure_options="$configure_options --with-boost=yes --with-boost-dir=$1/BOOST "
  configure_trilinos_options="$configure_trilinos_options -D TPL_ENABLE_Boost=ON  -D TPL_ENABLE_BoostLib=ON  -D Boost_INCLUDE_DIRS=$1/BOOST/include -D BoostLib_LIBRARY_DIRS=$1/BOOST/lib -D BoostLib_INCLUDE_DIRS=$1/BOOST/include"
fi

if [ -d "$1/MPI" ]; then
  configure_trilinos_options="$configure_trilinos_options -D MPI_BASE_DIR=$1/MPI "
  mpi_dir="$1/MPI"
else
  mpi_dir=$(dirname "$(dirname "$(which mpic++)")")
fi

### It seem that the PETSC --download-packege option has several problems and cannot produce
### a valid compilation command for most of the packages + it seem also that some library
### are compiled without optimization enabled, so we provide manual installation for that packages

if [ ! -d "$1/OPENBLAS" ]; then
  ./script/install_OPENBLAS.sh $1
  if [ $? -eq 0 ]; then
    configure_options="$configure_options --with-blas-lib=$1/OPENBLAS/lib/libopenblas.a --with-lapack-lib=$1/OPENBLAS/lib/libopenblas.a"
    configure_trilinos_options="$configure_trilinos_options -D TPL_ENABLE_BLAS=ON -D BLAS_LIBRARY_NAMES=openblas -D BLAS_LIBRARY_DIRS=$1/OPENBLAS/lib -D TPL_ENABLE_LAPACK=ON -D LAPACK_LIBRARY_NAMES=openblas -D LAPACK_LIBRARY_DIRS=$1/OPENBLAS/lib -D TPL_ENABLE_Netcdf=OFF -DTPL_ENABLE_GLM=OFF -D TPL_ENABLE_X11=OFF  "
    configure_options_superlu="$configure_options_superlu -Denable_blaslib=OFF  -DTPL_BLAS_LIBRARIES=$1/OPENBLAS/lib/libopenblas.a "
    configure_options_hypre="--with-blas-libs=-lopenblas --with-blas-lib-dirs=$1/OPENBLAS/lib --with-lapack-libs=-lopenblas  --with-lapack-lib-dirs=$1/OPENBLAS/lib "
  fi
else
    configure_options="$configure_options --with-blas-lib=$1/OPENBLAS/lib/libopenblas.a --with-lapack-lib=$1/OPENBLAS/lib/libopenblas.a"
    configure_trilinos_options="$configure_trilinos_options -D TPL_ENABLE_BLAS=ON -D BLAS_LIBRARY_NAMES=openblas -D BLAS_LIBRARY_DIRS=$1/OPENBLAS/lib -D TPL_ENABLE_LAPACK=ON -D LAPACK_LIBRARY_NAMES=openblas -D LAPACK_LIBRARY_DIRS=$1/OPENBLAS/lib -D TPL_ENABLE_Netcdf=OFF -DTPL_ENABLE_GLM=OFF -D TPL_ENABLE_X11=OFF  "
    configure_options_superlu="$configure_options_superlu -Denable_blaslib=OFF  -DTPL_BLAS_LIBRARIES=$1/OPENBLAS/lib/libopenblas.a "
    configure_options_hypre="--with-blas-libs=-lopenblas --with-blas-lib-dirs=$1/OPENBLAS/lib --with-lapack-libs=-lopenblas  --with-lapack-lib-dirs=$1/OPENBLAS/lib "
fi

if [ ! -d "$1/SUITESPARSE" ]; then
  CXX="$CXX" CC="$CC" FC="$FC" F77="$F77" ./script/install_SUITESPARSE.sh $1 $2
  if [ $? -eq 0 ]; then
    configure_options="$configure_options --with-suitesparse=yes --with-suitesparse-dir=$1/SUITESPARSE "
  fi
else
  configure_options="$configure_options --with-suitesparse=yes --with-suitesparse-dir=$1/SUITESPARSE "
fi

configure_options="$configure_options --download-scalapack --download-mumps"
configure_options="$configure_options --download-superlu_dist"
configure_options="$configure_options --download-hypre"

rm petsc-lite-3.10.2.tar.gz
rm -rf petsc-3.10.2
wget http://ppmcore.mpi-cbg.de/upload/petsc-lite-3.10.2.tar.gz
if [ $? -ne 0 ]; then
  echo -e "\033[91;5;1m FAILED! Installation requires an Internet connection \033[0m"
  exit 1
fi
tar -xf petsc-lite-3.10.2.tar.gz
cd petsc-3.10.2

if [ x"$CXX" != x"icpc" ]; then

  echo "./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc --with-cxx-dialect=C++11 $petsc_openmp  --with-mpi-dir=$mpi_dir $configure_options  --prefix=$1/PETSC --with-debugging=0"

  function haveProg() {
      [ -x "$(command -v $1)" ]
  }

  $python_command ./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc  --with-cxx-dialect=C++11 $petsc_openmp --with-mpi-dir=$mpi_dir $configure_options --prefix=$1/PETSC --with-debugging=0

else

  echo "./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc --with-cxx-dialect=C++11 $petsc_openmp  --with-mpi-dir=$mpi_dir $configure_options  --prefix=$1/PETSC --with-debugging=0"

  function haveProg() {
      [ -x "$(command -v $1)" ]
  }

  $python_command ./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc  --with-cxx-dialect=C++11 $petsc_openmp --with-mpi-dir=$mpi_dir $configure_options --prefix=$1/PETSC --with-debugging=0

fi

make all test
make install

# if empty remove the folder
if [ ! "$(ls -A $1/PETSC)" ]; then
   rm -rf $1/PETSC
else
   #Mark the installation
   echo 2 > $1/PETSC/version
   exit 0
fi

