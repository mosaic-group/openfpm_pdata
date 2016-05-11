#! /bin/bash

# check if the directory $1/EIGEN exist

if [ -d "$1/EIGEN" ]; then
  echo "EIGEN already installed"
  exit 0
fi

./script/install_OPENBLAS.sh $1 $2
if [ ! -d "$1/OPENBLAS"  ]; then
  exit 1
fi

./script/install_SUITESPARSE.sh $1 $2
if [ ! -d "$1/SUITESPARSE"  ]; then
  exit 1
fi

wget http://ppmcore.mpi-cbg.de/upload/eigen-3.2.7.tar.bz2
rm -rf eigen-eigen-b30b87236a1b
tar -xf eigen-3.2.7.tar.bz2

cd eigen-eigen-b30b87236a1b
mkdir $1/EIGEN/
mv Eigen $1/EIGEN/Eigen

cd ..
rm -rf eigen-eigen-b30b87236a1b

