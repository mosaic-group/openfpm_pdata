#! /bin/bash

source script/detect_gcc
source script/discover_os

discover_os

if [ x"$platform" == x"msys" ]; then
	pacman -Syu mingw-w64-x86_64-suitesparse
else

	# check if the directory $1/SUITESPARSE exist
	rm -rf SuiteSparse-5.7.2
	if [ -d "$1/SUITESPARSE"  -a -f "$1/SUITESPARSE/include/umfpack.h" ]; then
  		echo "SUITESPARSE is already installed"
  		exit 0
	fi

	rm SuiteSparse-5.7.2.tar.gz
	wget http://ppmcore.mpi-cbg.de/upload/SuiteSparse-5.7.2.tar.gz
	rm -rf SuiteSparse
	tar -xf SuiteSparse-5.7.2.tar.gz
	if [ $? != 0 ]; then
  		echo "Failed to download SuiteSparse"
  		exit 1
	fi
	cd SuiteSparse-5.7.2

	if [ x"$CXX" == x"icpc" ]; then
    		export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/$1/OPENBLAS/lib"
    		STS_LIB="-shared-intel -lrt -lifcore"
	fi

	export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$1/OPENBLAS/lib"

	if [ x"$platform" == x"cygwin" ]; then
    		export PATH="$PATH:$(pwd)/lib"
    		echo "$PATH"
	fi

	echo "Compiling SuiteSparse without CUDA (old variable $CUDA)"
	LDLIBS="$STS_LIB -lm" LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$1/OPENBLAS/lib"  make library -j $2 "CC=$CC" "CXX=$CXX" "CUDA=no" "BLAS=-L$1/OPENBLAS/lib -lopenblas -pthread" "LAPACK=-lopenblas"
	if [ $? != 0 ]; then
  		echo "Failed to compile SuiteSparse"
  		exit 1
	fi
	echo "Making library"
	make library "CC=$CC" "CXX=$CXX" "CUDA=no" "INSTALL=$1/SUITESPARSE" "INSTALL_LIB=$1/SUITESPARSE/lib" "INSTALL_INCLUDE=$1/SUITESPARSE/include" "BLAS=-L$1/OPENBLAS/lib -lopenblas -pthread" "LAPACK="
	echo "Making install"
	make install "CC=$CC" "CXX=$CXX" "CUDA=no" "INSTALL=$1/SUITESPARSE" "INSTALL_LIB=$1/SUITESPARSE/lib" "INSTALL_INCLUDE=$1/SUITESPARSE/include" "BLAS=-L$1/OPENBLAS/lib -lopenblas -pthread" "LAPACK="
	# Mark the installation

	rm SuiteSparse-5.7.2.tar.gz
	echo 2 > $1/SUITESPARSE/version
fi


