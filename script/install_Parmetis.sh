#! /bin/bash

source script/discover_os
discover_os

# check if the directory $1/PARMETIS exist

if [ -f "$1/PARMETIS/include/parmetis.h" ]; then
  echo "PARMETIS is already installed"
  exit 0
fi

## Remove old download
rm -rf parmetis-4.0.3
rm parmetis-4.0.3.tar.gz

#wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
wget http://ppmcore.mpi-cbg.de/upload/parmetis-4.0.3.tar.gz
tar -xf parmetis-4.0.3.tar.gz
cd parmetis-4.0.3
# Change to 64 bit parmetis

if [ x"$platform" == x"osx" ]; then
  sed_command=gsed
else
  sed_command=sed
fi

$sed_command -i "/#define\sIDXTYPEWIDTH\s32/c\#define IDXTYPEWIDTH 64" metis/include/metis.h

CFLAGS=-fPIC make config prefix=$1/PARMETIS
make -j $2
if [ $? -ne 0 ]; then
    echo "PARMETIS error installing"
    exit 0
fi
mkdir $1/PARMETIS
make install

#### Apply patch if we are on cygwin

if [ x"$platform" == x"cygwin" ]; then
  cd $1/PARMETIS/include
  wget http://openfpm.mpi-cbg.de/upload/parmetis_patch
  patch < parmetis_patch
fi

# Mark the installation
echo 3 > $1/PARMETIS/version

