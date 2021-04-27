#!/bin/bash

# check if the directory $1/MPI exist

if [ -d "$1/MPI" ]; then
        version=$(cat $1/MPI/version)
        if [ x"$version" != x"10"  ]; then
            echo -e "\033[1;34;5m  -------------------------------------------------------------------------------------- \033[0m"
            echo -e "\033[1;34;5m  MPICH has been updated to version 3.3.0, the component will be updated automatically      \033[0m"
            echo -e "\033[1;34;5m  -------------------------------------------------------------------------------------- \033[0m"
            sleep 5
            rm -rf $1/MPI/include
            rm -rf $1/MPI/lib
            rm -rf $1/MPI/bin
            rm -rf $1/MPI/etc
            rm -rf $1/MPI/share
            rm -rf $1/MPI
            rm -rf $1/HDF5
            rm -rf $1/ZLIB
            rm -rf $1/PARMETIS
            rm -rf $1/PETSC
            rm -rf $1/TRILINOS
            rm -rf $1/HYPRE
            rm -rf $1/MUMPS
            rm -rf $1/SUPERLU_DIST
        else
                echo "MPI already installed"
                exit 0
        fi
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
echo 10 > $1/MPI/version

