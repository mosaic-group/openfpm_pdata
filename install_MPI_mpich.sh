#!/bin/bash

# check if the directory $1/MPI exist

if [ -d "$1/MPI" ]; then
  echo "MPI already installed"
  exit 0
fi
rm -rf mpich-3.3
rm mpich-3.3.tar.gz
wget http://ppmcore.mpi-cbg.de/upload/mpich-3.3.tar.gz
tar -xvf mpich-3.3.tar.gz
cd mpich-3.3

#
#                  --disable-mca-dso \
#                 --disable-sysv-shmem \
#                 --enable-cxx-exceptions \
#                 --with-threads=posix \
#                 --without-cs-fs \
#                 --with-mpi-param_check=always \
#                 --enable-contrib-no-build=vt,libompitrace \
#
#--enable-mca-no-build=paffinity,installdirs-windows,timer-windows,shmem-sysv
#
#


if [ x"$6" == x"1" ]; then
   echo "Installing MPI with GPU support"
  ./configure --prefix=$1/MPI --enable-fortran CC=$3 CXX=$4 F77=$5 FC=$5
else
  echo "Installing MPI without GPU support"
  ./configure --prefix=$1/MPI --enable-fortran CC=$3 CXX=$4 F77=$5 FC=$5
fi
make -j $2
make install

# Mark the installation
echo 4 > $1/MPI/version

