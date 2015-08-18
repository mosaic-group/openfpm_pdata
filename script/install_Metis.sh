#! /bin/bash

# check if the directory ${HOME}/METIS exist

if [ -d "${HOME}/METIS" ]; then
  echo "METIS already installed"
  exit 0
fi


mkdir ${HOME}/METIS
wget http://ppmcore.mpi-cbg.de/upload/metis-5.1.0.tar.gz
tar -xf metis-5.1.0.tar.gz
cd metis-5.1.0
make config shared=1 prefix=${HOME}/METIS
make
make install

