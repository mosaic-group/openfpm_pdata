#!/bin/bash 

# check if the directory ${HOME}/MPI exist

if [ -d "${HOME}/MPI" ]; then
  echo "MPI already installed"
  exit 0
fi

mkdir ${HOME}/MPI
wget http://www.open-mpi.de/software/ompi/v1.8/downloads/openmpi-1.8.7.tar.bz2
tar -xvf openmpi-1.8.7.tar.bz2
cd openmpi-1.8.7
sh ./configure --prefix=${HOME}/MPI --enable-opal-multi-threads --enable-mpi-f90
make -j 4
make install
