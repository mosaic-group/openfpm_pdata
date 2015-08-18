#!/bin/bash 

# check if the directory ${HOME}/MPI exist

if [ -d "${HOME}/BOOST" ]; then
  echo "BOOST already installed"
  exit 0
fi

mkdir ${HOME}/BOOST
wget http://ppmcore.mpi-cbg.de/upload/boost_1_58_0.tar.bz2
tar -xvf boost_1_58_0.tar.bz2
cd boost_1_58_0
./bootstrap.sh
./b2 -j 4 install --prefix=$HOME/BOOST

