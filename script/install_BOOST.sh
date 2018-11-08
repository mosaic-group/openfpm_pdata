#!/bin/bash 

# check if the directory $1/BOOST exist

if [ -d "$1/BOOST" ]; then
  echo "BOOST already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/boost_1_68_0.tar.bz2
tar -xvf boost_1_68_0.tar.bz2
cd boost_1_68_0
./bootstrap.sh --with-toolset=$3
mkdir $1/BOOST
./b2 -j $2 install --prefix=$1/BOOST
rm -rf boost_1_68_0

