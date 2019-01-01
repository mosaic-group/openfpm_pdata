#!/bin/bash

# check if the directory $1/MPI exist

./script/download_MPI.sh

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


if [ x"$3" == x"1" ]; then
   echo "Installing MPI with GPU support"
  ./configure --with-cuda --prefix=$1/MPI --enable-mpi-fortran=yes CC=$4 CXX=$5 F77=$6 FC=$7
else
  echo "Installing MPI without GPU support"
  ./configure --prefix=$1/MPI --enable-mpi-fortran=yes CC=$4 CXX=$5 F77=$6 FC=$7
fi
make -j $2
make install

# Mark the installation
echo 4 > $1/MPI/version

