#! /bin/bash

wget http://ppmcore.mpi-cbg.de/upload/eigen-3.4.0.tar.bz2
tar -xf eigen-3.4.0.tar.bz2

cd eigen-3.4.0/
mkdir $1/EIGEN/
mv Eigen $1/EIGEN/Eigen
