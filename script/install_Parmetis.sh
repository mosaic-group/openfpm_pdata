#! /bin/bash

# check if the directory $1/PARMETIS exist

if [ -d "$1/PARMETIS" ]; then
  echo "PARMETIS already installed"
  exit 0
fi

## Remove old download
rm -rf parmetis-4.0.3

wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
tar -xf parmetis-4.0.3.tar.gz
cd parmetis-4.0.3
# Change to 64 bit parmetis

if [ x"$platform" == x"osx" ]; then
  sed_command=gsed
else
  sed_command=sed
fi

$sed_command -i "/#define\sIDXTYPEWIDTH\s32/c\#define IDXTYPEWIDTH 64" metis/include/metis.h

sed #define IDXTYPEWIDTH 32

make config prefix=$1/PARMETIS
make -j $2
if [ $? -ne 0 ]; then
    echo "PARMETIS error installing"
    exit 0
fi
mkdir $1/PARMETIS
make install

# Mark the installation
echo 1 > $1/PARMETIS/version

