#! /bin/bash

# check if the directory $1/METIS exist

if [ -d "$1/METIS" ]; then
  echo "METIS already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/metis-5.1.0.tar.gz
tar -xf metis-5.1.0.tar.gz
cd metis-5.1.0
make config shared=1 prefix=$1/METIS
make -j 4
mkdir $1/METIS
make install

