#!/bin/bash 

# check if the directory $1/MPI exist

if [ -d "$1/BOOST" ]; then
  echo "BOOST already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/boost_1_58_0.tar.bz2
tar -xvf boost_1_60_0.tar.bz2
cd boost_1_60_0
./bootstrap.sh
mkdir $1/BOOST
./b2 -j 4 install --prefix=$1/BOOST

