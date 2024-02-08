#! /bin/bash

rm -rf SuiteSparse
wget http://ppmcore.mpi-cbg.de/upload/SuiteSparse-5.7.2.tar.gz -O SuiteSparse-5.7.2.tar.gz
tar -xf SuiteSparse-5.7.2.tar.gz
cd SuiteSparse-5.7.2

if [ x"$CXX" == x"icpc" ]; then
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/$1/OPENBLAS/lib"
    STS_LIB="-shared-intel -lrt -lifcore"
fi

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$1/OPENBLAS/lib"

if [[ "$OSTYPE" == "cygwin" ]]; then
    export PATH="$PATH:$(pwd)/lib"
    echo "$PATH"
fi

echo "Compiling SuiteSparse without CUDA (old variable $CUDA)"
LDLIBS="$STS_LIB -lm" LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$1/OPENBLAS/lib"  make library "CC=$CC" "CXX=$CXX" "CUDA=no" "BLAS=-L$1/OPENBLAS/lib -lopenblas -pthread" "LAPACK=-lopenblas" -j $2

echo "Making library"
make library "CC=$CC" "CXX=$CXX" "CUDA=no" "INSTALL=$1/SUITESPARSE" "INSTALL_LIB=$1/SUITESPARSE/lib" "INSTALL_INCLUDE=$1/SUITESPARSE/include" "BLAS=-L$1/OPENBLAS/lib -lopenblas -pthread" "LAPACK=" -j $2
echo "Making install"
make install "CC=$CC" "CXX=$CXX" "CUDA=no" "INSTALL=$1/SUITESPARSE" "INSTALL_LIB=$1/SUITESPARSE/lib" "INSTALL_INCLUDE=$1/SUITESPARSE/include" "BLAS=-L$1/OPENBLAS/lib -lopenblas -pthread" "LAPACK="
