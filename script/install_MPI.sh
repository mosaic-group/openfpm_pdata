#!/bin/bash

# check if the directory $1/MPI exist

if [ -d "$1/MPI" ]; then
  echo "MPI already installed"
  exit 0
fi

./script/download_MPI.sh
cd openmpi-4.1.1

if [ -f ../mpi_add_options ]; then
	mpi_options=$(cat ../mpi_add_options)
	echo "Adding MPI options: $mpi_options"
fi

if [ x"$3" == x"1" ]; then
   echo "Installing MPI with GPU support"

   # Detect where is nvcc
   cuda_location=$(dirname $(dirname $(which nvcc)) )

   ./configure --with-hwloc=internal --with-libevent=internal $mpi_options --with-cuda=$cuda_location --prefix=$1/MPI --enable-mpi-fortran=yes CC=$4 CXX=$5 F77=$6 FC=$7 $8
else
   echo "Installing MPI without GPU support"
   ./configure --with-hwloc=internal --with-libevent=internal $mpi_options --prefix=$1/MPI --enable-mpi-fortran=yes CC=$4 CXX=$5 F77=$6 FC=$7 $8
fi
make -j $2
make install

# Mark the installation
echo 10 > $1/MPI/version

