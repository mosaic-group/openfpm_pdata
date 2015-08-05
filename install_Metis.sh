#! /bin/bash

mkdir ${HOME}/METIS
wget http://ppmcore.mpi-cbg.de/upload/metis-5.1.0.tar.gz
tar -xf metis-5.1.0.tar.gz
cd metis-5.1.0.tar.gz
sh ./configure --prefix={$HOME}/METIS
make
make install
