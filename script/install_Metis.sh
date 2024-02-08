#! /bin/bash

rm -rf metis-5.1.0
wget http://ppmcore.mpi-cbg.de/upload/metis-5.1.0.tar.gz -O metis-5.1.0.tar.gz
tar -xf metis-5.1.0.tar.gz
cd metis-5.1.0

if [[ "$OSTYPE" == "cygwin" ]]; then
	shared_opt="-DSHARED=OFF"
else
	shared_opt="-DSHARED=ON"
fi

cputype=$(uname -m | sed "s/\\ /_/g")
systype=$(uname -s)
BUILDDIR=build/$systype-$cputype
mkdir -p $BUILDDIR
cd $BUILDDIR
if [ "$#" -eq 4  ]; then
  echo "cmake ../../. $shared_opt -DGKLIB_PATH=../../GKlib -DCMAKE_INSTALL_PREFIX=$1/METIS -DCMAKE_C_COMPILER=$3 -DCMAKE_CXX_COMPILER=$4"
  cmake ../../. $shared_opt -DGKLIB_PATH=../../GKlib  -DCMAKE_INSTALL_PREFIX=$1/METIS -DCMAKE_C_COMPILER=$3 -DCMAKE_CXX_COMPILER=$4
else
  echo "cmake ../../. $shared_opt -DGKLIB_PATH=../../GKlib -DCMAKE_INSTALL_PREFIX=$1/METIS"
  cmake ../../. $shared_opt -DGKLIB_PATH=../../GKlib -DCMAKE_INSTALL_PREFIX=$1/METIS
fi
make -j $2
make install

# Mark the installation
echo 3 > $1/METIS/version

