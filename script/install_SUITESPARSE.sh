#! /bin/bash

source script/detect_gcc
source script/discover_os

discover_os

# check if the directory $1/SUITESPARSE exist

if [ -d "$1/SUITESPARSE" ]; then
  echo "SUITESPARSE already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/SuiteSparse-4.5.5.tar.gz
rm -rf SuiteSparse
tar -xf SuiteSparse-4.5.5.tar.gz
if [ $? != 0 ]; then
  echo "Fail to download SuiteSparse"
  exit 1
fi
cd SuiteSparse

# configuration

#if [ x"$platform" = x"osx"  ]; then
#    # installation for OSX

#    sed -i "" -e "s| LAPACK = -llapack|LAPACK = |" SuiteSparse_config/SuiteSparse_config_Mac.mk
#    sed -i "" -e "s| BLAS = -lopenblas|BLAS = -L"$1"/OPENBLAS/lib -lopenblas|" SuiteSparse_config/SuiteSparse_config_Mac.mk

    ### Overwrite SuiteSparse_config.mk

#    rm SuiteSparse_config/SuiteSparse_config.mk
#    mv SuiteSparse_config/SuiteSparse_config_Mac.mk SuiteSparse_config/SuiteSparse_config.mk

#else
    
    # Installation for linux

    if [ x"$CXX" == x"icpc" ]; then
      export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/$1/OPENBLAS/lib"
      STS_LIB="-shared-intel -lrt -lifcore"
    fi
#      sed -i "/\sLIB\s=\s-lm\s-lrt/c\LIB = -shared-intel -lm -lrt -lifcore" SuiteSparse_config/SuiteSparse_config.mk
#      sed -i "/\sLAPACK\s=\s-llapack/c\LAPACK = " SuiteSparse_config/SuiteSparse_config.mk
#      sed -i "/\sBLAS\s=\s\-lopenblas/c\BLAS = -L$1/OPENBLAS/lib -lopenblas -lpthread" SuiteSparse_config/SuiteSparse_config.mk
#    else
#      sed -i "/\sLAPACK\s=\s-llapack/c\LAPACK = " SuiteSparse_config/SuiteSparse_config.mk
#      sed -i "/\sBLAS\s=\s\-lopenblas/c\BLAS = -L$1/OPENBLAS/lib -lopenblas -lpthread" SuiteSparse_config/SuiteSparse_config.mk
#    fi

#fi

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$1/OPENBLAS/lib"

echo "Compiling SuiteSparse without CUDA (old variable $CUDA)"
make "CUDA=no" "BLAS=-L$1/OPENBLAS/lib -lopenblas" "LAPACK="
if [ $? != 0 ]; then
  echo "Fail to compile SuiteSparse"
  exit 1
fi
make install "CUDA=no" "INSTALL=$1/SUITESPARSE" "INSTALL_LIB=$1/SUITESPARSE/lib" "INSTALL_INCLUDE=$1/SUITESPARSE/include" "BLAS=-L$1/OPENBLAS/lib -lopenblas" "LAPACK="
# Mark the installation
echo 1 > $1/SUITESPARSE/version
rm -rf SuiteSparse
rm SuiteSparse-4.5.5.tar.gz
