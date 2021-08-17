#!/bin/bash

if [ -d "$1/CMAKE" ]; then
   exit 0
fi

wget https://github.com/Kitware/CMake/releases/download/v3.20.3/cmake-3.20.3.tar.gz
tar -xvf cmake-3.20.3.tar.gz
cd cmake-3.20.3

./bootstrap --prefix=$1/CMAKE
make
make install


