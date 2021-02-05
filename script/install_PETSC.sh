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

function test_configure_options() {
  cd petsc-3.14.2
  $python_command ./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc  --with-cxx-dialect=C++11 $petsc_openmp --with-mpi-dir=$mpi_dir $configure_options2 --with-debugging=0
  error=$?
  cd ..
}

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

#### Download and uncompress petsc

rm petsc-lite-3.14.2.tar.gz
rm -rf petsc-3.14.2
wget http://ppmcore.mpi-cbg.de/upload/petsc-lite-3.14.2.tar.gz
if [ $? -ne 0 ]; then
  echo -e "\033[91;5;1m FAILED! Installation requires an Internet connection \033[0m"
  exit 1
fi
tar -xf petsc-lite-3.14.2.tar.gz

####

## If some dependencies has been installed feed them to PETSC

MUMPS_extra_libs=""

configure_options=""


configure_options="$configure_options --download-metis --download-parmetis"

if [ -d "$1/BOOST" ]; then

  ### We check incrementaly the options
  configure_options="$configure_options --with-boost=yes --with-boost-dir=$1/BOOST "
fi

if [ -d "$1/MPI" ]; then
  mpi_dir="$1/MPI"
else
  mpi_dir=$(dirname "$(dirname "$(which mpic++)")")
fi


if [ ! -d "$1/OPENBLAS" ]; then
  CXX="$CXX" CC="$CC" FC="$FC" F77="$F77" ./script/install_OPENBLAS.sh $1
  if [ $? -eq 0 ]; then
    configure_options="$configure_options --with-blas-lib=$1/OPENBLAS/lib/libopenblas.a --with-lapack-lib=$1/OPENBLAS/lib/libopenblas.a"
  fi
else
    configure_options="$configure_options --with-blas-lib=$1/OPENBLAS/lib/libopenblas.a --with-lapack-lib=$1/OPENBLAS/lib/libopenblas.a"
fi

if [ -d "$1/SUITESPARSE"  -a -f "$1/SUITESPARSE/include/umfpack.h" ]; then
  CXX="$CXX" CC="$CC" FC="$FC" F77="$F77" ./script/install_SUITESPARSE.sh $1 $2
fi

#### OK here we check if we can configure work with SUITESPARSE
if [ -d "$1/SUITESPARSE"  -a -f "$1/SUITESPARSE/include/umfpack.h" ]; then
  echo "Testing if PETSC work with SUITESPARSE"
  configure_options2="$configure_options --with-suitesparse=yes --with-suitesparse-dir=$1/SUITESPARSE "
  test_configure_options

  if [ $error -eq 0 ]; then
    echo "SUITESPARSE work with PETSC"
    configure_options="$configure_options --with-suitesparse=yes --with-suitesparse-dir=$1/SUITESPARSE "
  fi
fi

configure_options2="$configure_options --download-scalapack"
test_configure_options

if [ $error -eq 0 ]; then
  echo "SCALAPACK work with PETSC"
  configure_options="$configure_options --download-scalapack "
fi

configure_options2="$configure_options --download-mumps"
test_configure_options

if [ $error -eq 0 ]; then
  echo "MUMPS work with PETSC"
  configure_options="$configure_options --download-mumps"
fi


#### OK here we check if we can configure work with SUITESPARSE
echo "Testing if PETSC work with SUPERLU"
configure_options2="$configure_options --download-superlu_dist "
test_configure_options

if [ $error -eq 0 ]; then
  echo "SUITESPARSE work with PETSC"
  configure_options="$configure_options --download-superlu_dist "
fi

configure_options2="$configure_options --download-hypre"
test_configure_options

if [ $error -eq 0 ]; then
  echo "HYPRE work with PETSC"
  configure_options="$configure_options --download-hypre"
fi


rm petsc-lite-3.14.2.tar.gz
rm -rf petsc-3.14.2
wget http://ppmcore.mpi-cbg.de/upload/petsc-lite-3.14.2.tar.gz
if [ $? -ne 0 ]; then
  echo -e "\033[91;5;1m FAILED! Installation requires an Internet connection \033[0m"
  exit 1
fi
tar -xf petsc-lite-3.14.2.tar.gz
cd petsc-3.14.2

if [ x"$CXX" != x"icpc" ]; then

  echo "./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc --with-cxx-dialect=C++11 $petsc_openmp  --with-mpi-dir=$mpi_dir $configure_options  --prefix=$1/PETSC --with-debugging=0"

  function haveProg() {
      [ -x "$(command -v $1)" ]
  }

  if [ x"$platform" != x"cygwin" ]; then
  	$python_command ./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc  --with-cxx-dialect=C++11 $petsc_openmp --with-mpi-dir=$mpi_dir $configure_options --prefix=$1/PETSC --with-debugging=0
  else
	echo "Sorry PETSC installation in not supported on CYGWIN"
  fi

else

  echo "./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc --with-cxx-dialect=C++11 $petsc_openmp  --with-mpi-dir=$mpi_dir $configure_options  --prefix=$1/PETSC --with-debugging=0"

  function haveProg() {
      [ -x "$(command -v $1)" ]
  }

  $python_command ./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc  --with-cxx-dialect=C++11 $petsc_openmp --with-mpi-dir=$mpi_dir $configure_options --prefix=$1/PETSC --with-debugging=0
fi

make all
make install

# if empty remove the folder
if [ ! "$(ls -A $1/PETSC)" ]; then
   rm -rf $1/PETSC
else
   #Mark the installation
   echo 5 > $1/PETSC/version
   exit 0
fi

