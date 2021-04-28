#! /bin/bash

source script/discover_os

discover_os

# check if the directory $1/HDF5 exist

if [ -d "$1/HDF5" -a -f "$1/HDF5/include/hdf5.h" ]; then
  echo "HDF5 is already installed"
  exit 0
fi

if [ ! -d "$1/ZLIB"  -a x"$platform" != x"cygwin" -a x"$platform" != x"msys" ]; then
  rm zlib-1.2.11.tar.gz
  rm -rf zlib-1.2.11
  wget https://zlib.net/zlib-1.2.11.tar.gz
  if [ $? -ne 0 ]; then
    echo -e "\033[91;5;1m FAILED! Installation requires an Internet connection \033[0m"
    exit 1
  fi
  tar -xf zlib-1.2.11.tar.gz
  cd zlib-1.2.11

  CC=mpicc ./configure --prefix=$1/ZLIB
  make -j $2
  cd ..
  if [ $? -eq 0 ]; then
    make check install
  else
    echo -e "\033[91;5;1m ZLIB Installation FAILED \033[0m"
    exit 1
  fi

else
  echo "ZLIB is already installed"
fi

if [ x"$platform" == x"msys" ]; then
	wget http://ppmcore.mpi-cbg.de/upload/hdf5-1.12.0_msys.tar.bz2
	tar -xf hdf5-1.12.0_msys.tar.bz2
	cd hdf5-1.12.0
else
	wget http://ppmcore.mpi-cbg.de/upload/hdf5-1.10.7.tar.gz
	tar -xf hdf5-1.10.7.tar.gz
	cd 1.10.7
fi


if [ x"$platform" == x"cygwin" ]; then
        CC=mpicc ./configure --enable-parallel --prefix=$1/HDF5
	make CFLAGS=-D_POSIX_C_SOURCE -j $2
elif [ x"$platform" == x"msys" ]; then
	mkdir build
	cd build
	cmake -Wno-dev     -G"MSYS Makefiles"  -DHDF5_BUILD_FORTRAN=OFF   -DBUILD_SHARED_LIBS=ON     -DBUILD_TESTING=OFF     -DCMAKE_BUILD_TYPE=Release     -DCMAKE_SKIP_RPATH=ON     -DHDF5_BUILD_HL_LIB=ON     -DHDF5_BUILD_CPP_LIB=OFF     -DHDF5_BUILD_FORTRAN=ON     -DHDF5_BUILD_TOOLS=ON     -DHDF5_ENABLE_DEPRECATED_SYMBOLS=ON     -DHDF5_ENABLE_SZIP_SUPPORT=ON     -DHDF5_ENABLE_Z_LIB_SUPPORT=ON     -DHDF5_INSTALL_CMAKE_DIR="lib/cmake" -DHDF5_ENABLE_PARALLEL=ON   -DHDF5_INSTALL_DATA_DIR="share/doc/hdf5" ../. -DCMAKE_INSTALL_PREFIX=$1/HDF5
	make -j $2
else
        CC=mpicc ./configure --with-zlib=$1/ZLIB --enable-parallel --prefix=$1/HDF5
	make -j $2
fi
mkdir $1/HDF5
make install
if [ $? -ne 0 ]; then
    echo "HDF5 error installing"
    exit 0
fi
echo 3 > $1/HDF5/version
