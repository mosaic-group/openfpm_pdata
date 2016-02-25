#!/bin/bash 

# check if the directory $1/MPI exist

if [ -d "$1/MPI" ]; then
  echo "MPI already installed"
  exit 0
fi

wget http://www.open-mpi.de/software/ompi/v1.8/downloads/openmpi-1.8.7.tar.bz2
tar -xvf openmpi-1.8.7.tar.bz2
cd openmpi-1.8.7

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

sh ./configure --prefix=$1/MPI --enable-opal-multi-threads --enable-mpi-f90 $2 $3
make -j 4
mkdir $1/MPI
make install
