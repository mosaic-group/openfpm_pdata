#! /bin/bash

# check if the directory $1/EIGEN exist

if [ -d "$1/EIGEN" ]; then
  echo "EIGEN is already installed"
  exit 0
fi

./script/install_OPENBLAS.sh $1 $2
if [ ! -d "$1/OPENBLAS"  ]; then
  exit 1
fi

CXX="$CXX" CC="$CC" FC="$FC" F77="$F77" ./script/install_SUITESPARSE.sh $1 $2
if [ ! -d "$1/SUITESPARSE"  ]; then
  exit 1
fi

rm -rf eigen-3.3.5.tar.bz2
wget http://ppmcore.mpi-cbg.de/upload/eigen-3.3.5.tar.bz2
rm -rf eigen-eigen-b3f3d4950030/
tar -xf eigen-3.3.5.tar.bz2

cd eigen-eigen-b3f3d4950030/
mkdir $1/EIGEN/
mv Eigen $1/EIGEN/Eigen

cd ..
rm -rf eigen-eigen-b3f3d4950030/

touch $1/EIGEN/signature_of_eigen3_matrix_library

# Mark the installation
echo 2 > $1/EIGEN/version
