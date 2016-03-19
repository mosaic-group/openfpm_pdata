#! /bin/bash

# check if the directory $1/PARMETIS exist

if [ -d "$1/PARMETIS" ]; then
  echo "PARMETIS already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/parmetis-4.0.3.tar.gz
tar -xf parmetis-4.0.3.tar.gz
cd parmetis-4.0.3
make config prefix=$1/PARMETIS
make -j 4
if [ $? -ne 0 ]; then
    echo "PARMETIS error installing"
    exit 0
fi
mkdir $1/PARMETIS
make install

# Mark the installation
echo 1 > $1/PARMETIS/version

