#!/bin/bash

# check if the directory $1/MPI exist

if [ -d "$1/MPI" ]; then
  echo "MPI already installed"
  exit 0
fi
rm -rf openmpi-3.1.3
rm openmpi-3.1.3.tar.gz
wget http://ppmcore.mpi-cbg.de/upload/openmpi-3.1.3.tar.gz
tar -xvf openmpi-3.1.3.tar.gz

