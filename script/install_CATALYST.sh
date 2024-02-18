#! /bin/bash

source script/discover_os

discover_os

# check if the directory $1/CATALYST exist

if [ -d "$1/CATALYST" -a -f "$1/CATALYST/include/catalyst-2.0/catalyst.h" ]; then
  echo "CATALYST is already installed"
  exit 0
fi

git clone https://gitlab.kitware.com/paraview/catalyst.git
cd catalyst

mkdir build && cd build
cmake ../. -DCMAKE_INSTALL_PREFIX=$1/CATALYST
make -j $2

make install
cd ../..
rm -rf catalyst
if [ $? -ne 0 ]; then
    echo "CATALYST error installing"
    exit 0
fi
echo 1 > $1/CATALYST/version