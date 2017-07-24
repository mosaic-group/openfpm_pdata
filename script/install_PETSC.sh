#! /bin/bash

# check if the directory $1/PETSC exist

CXX=$4
CC=$3
F77=$5
FC=$6

if [ -d "$1/PETSC" ]; then
  echo "PETSC is already installed"
  exit 0
fi

# Detect gcc pr clang

source script/discover_os
discover_os

##### if we are on osx we use gsed

ldflags_petsc=
if [ x"$platform" == x"osx" ]; then
  sed_command=gsed
  ldflags_petsc=
else
  sed_command=sed
  ldflags_petsc=
fi


####

## If some dependencies has been installed feed them to PETSC

MUMPS_extra_libs=""

configure_options=""
configure_options_scalapack=""
configure_options_superlu=""
configure_trilinos_options=" -D TPL_ENABLE_MPI=ON "
configure_options_hypre=""

### Here we install OpenBLAS and SUITESPARSE


if [ -d "$1/PARMETIS" ]; then
  configure_options="$configure_options --with-parmetis=yes  --with-parmetis-dir=$1/PARMETIS "
  configure_options_superlu="-DTPL_PARMETIS_INCLUDE_DIRS=$1/PARMETIS/include;$1/METIS/include -DTPL_PARMETIS_LIBRARIES=$1/PARMETIS/lib/libparmetis.a;$1/METIS/lib/libmetis.so $configure_options_superlu"
fi

if [ -d "$1/METIS" ]; then
  configure_options="$configure_options --with-metis=yes --with-metis-dir=$1/METIS  "
fi

if [ -d "$1/BOOST" ]; then
  configure_options="$configure_options --with-boost=yes --with-boost-dir=$1/BOOST "
  configure_trilinos_options="$configure_trilinos_options -D TPL_ENABLE_Boost=ON  -D TPL_ENABLE_BoostLib=ON  -D Boost_INCLUDE_DIRS=$1/BOOST/include -D BoostLib_LIBRARY_DIRS=$1/BOOST/lib -D BoostLib_INCLUDE_DIRS=$1/BOOST/include"
fi

if [ -d "$1/MPI" ]; then
  configure_trilinos_options="$configure_trilinos_options -D MPI_BASE_DIR=$1/MPI "
  mpi_dir="$1/MPI"
else
  mpi_dir=$(dirname "$(dirname "$(which mpic++)")")
fi

### It seem that the PETSC --download-packege option has several problems and cannot produce
### a valid compilation command for most of the packages + it seem also that some library
### are compiled without optimization enabled, so we provide manual installation for that packages

if [ ! -d "$1/OPENBLAS" ]; then
  ./script/install_OPENBLAS.sh $1
  if [ $? -eq 0 ]; then
    configure_options="$configure_options --with-blas-lib=$1/OPENBLAS/lib/libopenblas.a --with-lapack-lib=$1/OPENBLAS/lib/libopenblas.a"
    configure_trilinos_options="$configure_trilinos_options -D TPL_ENABLE_BLAS=ON -D BLAS_LIBRARY_NAMES=openblas -D BLAS_LIBRARY_DIRS=$1/OPENBLAS/lib -D TPL_ENABLE_LAPACK=ON -D LAPACK_LIBRARY_NAMES=openblas -D LAPACK_LIBRARY_DIRS=$1/OPENBLAS/lib -D TPL_ENABLE_Netcdf=OFF -DTPL_ENABLE_GLM=OFF -D TPL_ENABLE_X11=OFF  "
    configure_options_superlu="$configure_options_superlu -Denable_blaslib=OFF  -DTPL_BLAS_LIBRARIES=$1/OPENBLAS/lib/libopenblas.a "
    configure_options_hypre="--with-blas-libs=-lopenblas --with-blas-lib-dirs=$1/OPENBLAS/lib --with-lapack-libs=-lopenblas  --with-lapack-lib-dirs=$1/OPENBLAS/lib "
    configure_options_scalapack="$configure_options_scalapack -D LAPACK_LIBRARIES=$1/OPENBLAS/lib/libopenblas.a -D BLAS_LIBRARIES=$1/OPENBLAS/lib/libopenblas.a"
  fi
else
    configure_options="$configure_options --with-blas-lib=$1/OPENBLAS/lib/libopenblas.a --with-lapack-lib=$1/OPENBLAS/lib/libopenblas.a"
    configure_trilinos_options="$configure_trilinos_options -D TPL_ENABLE_BLAS=ON -D BLAS_LIBRARY_NAMES=openblas -D BLAS_LIBRARY_DIRS=$1/OPENBLAS/lib -D TPL_ENABLE_LAPACK=ON -D LAPACK_LIBRARY_NAMES=openblas -D LAPACK_LIBRARY_DIRS=$1/OPENBLAS/lib -D TPL_ENABLE_Netcdf=OFF -DTPL_ENABLE_GLM=OFF -D TPL_ENABLE_X11=OFF  "
    configure_options_superlu="$configure_options_superlu -Denable_blaslib=OFF  -DTPL_BLAS_LIBRARIES=$1/OPENBLAS/lib/libopenblas.a "
    configure_options_hypre="--with-blas-libs=-lopenblas --with-blas-lib-dirs=$1/OPENBLAS/lib --with-lapack-libs=-lopenblas  --with-lapack-lib-dirs=$1/OPENBLAS/lib "
    configure_options_scalapack="$configure_options_scalapack -D LAPACK_LIBRARIES=$1/OPENBLAS/lib/libopenblas.a -D BLAS_LIBRARIES=$1/OPENBLAS/lib/libopenblas.a"
fi

if [ ! -d "$1/SUITESPARSE" ]; then
  ./script/install_SUITESPARSE.sh $1
  if [ $? -eq 0 ]; then
    configure_options="$configure_options --with-suitesparse=yes --with-suitesparse-dir=$1/SUITESPARSE "
  fi
else
  configure_options="$configure_options --with-suitesparse=yes --with-suitesparse-dir=$1/SUITESPARSE "
fi

# Install NETCFD
if [ -d "$1/NETCDF" ]; then
  echo "NETCDF is already installed"
  configure_options="$configure_options --with-netcdf=yes -with-netcdf-dir=$1/NETCDF --with-hdf5=yes --with-hdf5-dir=$1/HDF5 "
else
  if [ -d "$1/HDF5" ]; then
     configure_options="$configure_options --with-hdf5=yes --with-hdf5-dir=$1/HDF5  "
  else
     ./script/install_HDF5.sh $1 $2
     if [ $? -eq 0 ]; then
        configure_options="$configure_options --with-hdf5=yes --with-hdf5-dir=$1/HDF5  "
     fi
  fi

  rm netcdf-4.4.1.1.tar.gz
  rm -rf netcdf-4.4.1.1
  wget http://ppmcore.mpi-cbg.de/upload/netcdf-4.4.1.1.tar.gz
  if [ $? -ne 0 ]; then
    echo -e "\033[91;5;1m FAILED! Installation requires an Internet connection \033[0m"
    exit 1
  fi
  tar -xf netcdf-4.4.1.1.tar.gz
  cd netcdf-4.4.1.1

  if [ x"$platform" == x"cygwin" ]; then
    ./configure CC=mpicc CPPFLAGS="-I$1/HDF5/include -I$1/ZLIB/include " LDFLAGS="-L$1/HDF5/lib -L$1/ZLIB/lib" --disable-netcdf-4 --disable-dap --disable-shared --prefix=$1/NETCDF
  else
    ./configure CC=mpicc CPPFLAGS="-I$1/HDF5/include -I$1/ZLIB/include " LDFLAGS="-L$1/HDF5/lib -L$1/ZLIB/lib" --disable-dap --disable-shared --prefix=$1/NETCDF
  fi
  make -j $2

  if [ $? -eq 0 ]; then
    make install
    configure_options="$configure_options --with-netcdf=yes -with-netcdf-dir=$1/NETCDF "
  else
    echo -e "\033[91;5;1m FAILED! NETCDF Installation \033[0m"
    exit 1
  fi
fi

if [ ! -d "$1/TRILINOS" ]; then
  rm trilinos-12.10.1-Source.tar.bz2
  rm -rf trilinos-12.10.1-Source
  wget http://ppmcore.mpi-cbg.de/upload/trilinos-12.10.1-Source.tar.bz2
  if [ $? -ne 0 ]; then
    echo -e "\033[91;5;1m FAILED! Installation requires an Internet connection \033[0m"
    exit 1
  fi
  tar -xf trilinos-12.10.1-Source.tar.bz2
  cd trilinos-12.10.1-Source
  mkdir build
  cd build

  ### On clang we have no openMP
  petsc_openmp=""
  if [ x"$CXX" == x"clang++" ]; then
    conf_trl_openmp="-D Trilinos_ENABLE_OpenMP=OFF"
  elif [ x"$CXX" == x"icpc" ]; then

    configure_trilinos_options="$configure_trilinos_options -D Trilinos_ENABLE_Xpetra=OFF -D Trilinos_ENABLE_Amesos2=OFF -D Trilinos_ENABLE_Ifpack2=OFF -D Trilinos_ENABLE_Teko=OFF "
  else
    conf_trl_openmp="-D Trilinos_ENABLE_OpenMP=ON"
#    petsc_openmp="--with-openmp=yes"
  fi

  cmake -D CMAKE_INSTALL_PREFIX:PATH=$1/TRILINOS -D CMAKE_BUILD_TYPE=RELEASE $conf_trl_openmp -D Trilinos_ENABLE_TESTS=OFF  -D Trilinos_ENABLE_ALL_PACKAGES=ON $configure_trilinos_options  ../.

  make -j $2
  if [ $? -eq 0 ]; then
    make install
    # Mark the installation
    echo 1 > $1/TRILINOS/version
    configure_options="$configure_options --with-trilinos=yes -with-trilinos-dir=$1/TRILINOS"
  fi
else
  echo "Trilinos is already installed"
  configure_options="$configure_options --with-trilinos=yes -with-trilinos-dir=$1/TRILINOS"
fi

### Scalapack installation

if [ ! -d "$1/SCALAPACK" ]; then
  rm scalapack-2.0.2.tgz
  rm -rf scalapack-2.0.2
  wget http://ppmcore.mpi-cbg.de/upload/scalapack-2.0.2.tgz
  if [ $? -ne 0 ]; then
    echo -e "\033[91;5;1m FAILED! Installation requires an Internet connection \033[0m"
    exit 1
  fi
  tar -xf scalapack-2.0.2.tgz
  cd scalapack-2.0.2
  mkdir build
  cd build
  echo "cmake -D CMAKE_EXE_LINKER_FLAGS=-pthread  -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_Fortran_FLAGS_RELEASE=-fpic -D MPI_C_COMPILER_FLAGS=-fpic -D MPI_Fortran_COMPILER_FLAGS=-fpic -D CMAKE_C_FLAGS=-fpic -D CMAKE_INSTALL_PREFIX="$1/SCALAPACK" $configure_options_scalapack ../."
  cmake -D CMAKE_EXE_LINKER_FLAGS=-pthread  -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_Fortran_FLAGS_RELEASE=-fpic -D MPI_C_COMPILER_FLAGS=-fpic -D MPI_Fortran_COMPILER_FLAGS=-fpic -D CMAKE_C_FLAGS=-fpic -D CMAKE_INSTALL_PREFIX="$1/SCALAPACK" $configure_options_scalapack ../.
  make -j $2
  if [ $? -eq 0 ]; then
    make install
    configure_options="$configure_options --with-scalapack=yes -with-scalapack-dir=$1/SCALAPACK"
  fi
else
  echo "Scalapack is already installed"
  configure_options="$configure_options --with-scalapack=yes -with-scalapack-dir=$1/SCALAPACK"
fi

### MUMPS installation
if [ x"$CXX" != x"icpc" ]; then
  if [ ! -d "$1/MUMPS" ]; then
    rm MUMPS_5.0.1.tar.gz
    rm -rf MUMPS_5.0.1
    wget http://ppmcore.mpi-cbg.de/upload/MUMPS_5.0.1.tar.gz
    if [ $? -ne 0 ]; then
      echo -e "\033[91;5;1m FAILED! Installation requires an Internet connection \033[0m"
      exit 1
    fi
    tar -xf MUMPS_5.0.1.tar.gz
    cd MUMPS_5.0.1
    cp Make.inc/Makefile.inc.generic Makefile.inc

    # Installation for linux

    $sed_command -i "/CC\s\+=\scc/c\CC = mpicc" Makefile.inc
    $sed_command -i "/FC\s\+=\sf90/c\FC = mpif90" Makefile.inc
    $sed_command -i "/FL\s\+=\sf90/c\FL = mpif90" Makefile.inc

    $sed_command -i "/SCALAP\s\+=\s-lscalapack\s-lblacs/c\SCALAP = -L$1/SCALAPACK/lib -L$1/OPENBLAS/lib -lscalapack" Makefile.inc
    $sed_command -i "/LIBBLAS\s\+=\s\-lopenblas/c\LIBBLAS = -lopenblas" Makefile.inc

    $sed_command -i "/OPTF\s\+=\s\-O/c\OPTF = -fpic -O3" Makefile.inc
    $sed_command -i "/OPTC\s\+=\s\-O\s-I./c\OPTC = -fpic -O3 -I." Makefile.inc
    $sed_command -i "/OPTL\s\+=\s\-O/c\OPTL = -fpic -O3" Makefile.inc

    $sed_command -i "/LIBBLAS\s=\s-lblas/c\LIBBLAS = -lopenblas" Makefile.inc

    $sed_command -i "/INCPAR\s\+=\s\-I\/usr\/include/c\INCPAR =" Makefile.inc
    $sed_command -i "/LIBPAR\s\+=\s\$(SCALAP)\s\-L\/usr\/lib\s\-lmpi/c\LIBPAR = \$(SCALAP)" Makefile.inc

    make -j $2
    if [ $? -eq 0 ]; then
      ## Copy LIB and include in the target directory

      mkdir $1/MUMPS
      cp -r include $1/MUMPS
      cp -r lib $1/MUMPS

      MUMPS_extra_lib="-L$1/MUMPS/lib -ldmumps -lmumps_common -lpord -pthread "
      configure_options="$configure_options --with-mumps=yes --with-mumps-include=$1/MUMPS/include"

    fi
  else
    echo "MUMPS is already installed"
    MUMPS_extra_lib="-L$1/MUMPS/lib -ldmumps -lmumps_common -lpord -pthread "
    configure_options="$configure_options --with-mumps=yes --with-mumps-include=$1/MUMPS/include"
  fi
fi

## SuperLU installation

if [ ! -d "$1/SUPERLU_DIST" ]; then
  rm superlu_dist_5.1.3.tar.gz
  rm -rf SuperLU_DIST_5.1.3
  wget http://ppmcore.mpi-cbg.de/upload/superlu_dist_5.1.3.tar.gz
  if [ $? -ne 0 ]; then
    echo -e "\033[91;5;1m FAILED! Installation requires an Internet connection \033[0m"
    exit 1
  fi
  tar -xf superlu_dist_5.1.3.tar.gz
  cd SuperLU_DIST_5.1.3

  mkdir build
  cd build

  if [ x"$platform" == x"cygwin" ]; then
    cmake .. -DCMAKE_C_FLAGS="-fPIC -std=c99 "  -DTPL_BLAS_LIBRARIES="$1/OPENBLAS/lib/libopenblas.a"  -DCMAKE_INSTALL_PREFIX="$1/SUPERLU_DIST"  -DTPL_PARMETIS_INCLUDE_DIRS="$1/PARMETIS/include/;$1/METIS/include/" -DTPL_PARMETIS_LIBRARIES="$1/PARMETIS/lib/libparmetis.a;$1/METIS/lib/libmetis.dll.a;-lmpi;-lopen-rte;-lopen-pal"
  else
    cmake .. -DCMAKE_C_FLAGS="-fPIC -std=c99 "  -DTPL_BLAS_LIBRARIES="$1/OPENBLAS/lib/libopenblas.a"  -DCMAKE_INSTALL_PREFIX="$1/SUPERLU_DIST"  -DTPL_PARMETIS_INCLUDE_DIRS="$1/PARMETIS/include/;$1/METIS/include/" -DTPL_PARMETIS_LIBRARIES="$1/PARMETIS/lib/libparmetis.a;$1/METIS/lib/libmetis.so"
  fi

  # Installation for linux

  make
  if [ $? -eq 0 ]; then
     make install
     echo 1 > $1/SUPERLU_DIST/version

    if [ x"$CXX" == x"icpc" ]; then
      configure_options="$configure_options"
    else
      configure_options="$configure_options --with-superlu_dist=yes --with-superlu_dist-lib=$1/SUPERLU_DIST/lib/libsuperlu_dist.a --with-superlu_dist-include=$1/SUPERLU_DIST/include/"
    fi
  fi

else
  echo "SUPERLU is already installed"
  if [ x"$CXX" == x"icpc" ]; then
    configure_options="$configure_options"
  else
    configure_options="$configure_options --with-superlu_dist=yes --with-superlu_dist-lib=$1/SUPERLU_DIST/lib/libsuperlu_dist.a --with-superlu_dist-include=$1/SUPERLU_DIST/include/"
  fi
fi

## HYPRE installation

if [ ! -d "$1/HYPRE" ]; then
  rm hypre-2.11.2.tar.gz
  rm -rf hypre-2.11.2
  wget http://ppmcore.mpi-cbg.de/upload/hypre-2.11.2.tar.gz
  if [ $? -ne 0 ]; then
    echo -e "\033[91;5;1m FAILED! Installation requires an Internet connection \033[0m"
    exit 1
  fi
  tar -xf hypre-2.11.2.tar.gz
  cd hypre-2.11.2

  cd src

  ./configure CC=mpicc CXX=mpic++ CFLAGS=-fpic  $configure_options_hypre --prefix=$1/HYPRE
  make -j $2
  if [ $? -eq 0 ]; then
    make install
    echo 1 > $1/HYPRE/version
    configure_options="$configure_options --with-hypre=yes -with-hypre-dir=$1/HYPRE"
  fi

else
  echo "HYPRE is already installed"
  configure_options="$configure_options --with-hypre=yes -with-hypre-dir=$1/HYPRE"
fi

rm petsc-lite-3.7.6.tar.gz
rm -rf petsc-3.7.6
wget http://ppmcore.mpi-cbg.de/upload/petsc-lite-3.7.6.tar.gz
if [ $? -ne 0 ]; then
  echo -e "\033[91;5;1m FAILED! Installation requires an Internet connection \033[0m"
  exit 1
fi
tar -xf petsc-lite-3.7.6.tar.gz
cd petsc-3.7.6

echo "./configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc --with-cxx-dialect=C++11 $petsc_openmp  --with-mpi-dir=$mpi_dir $configure_options --with-mumps-lib="$MUMPS_extra_lib"  --prefix=$1/PETSC --with-debugging=0"

function haveProg() {
    [ -x "$(command -v $1)" ]
}

if haveProg python2; then
  python_command=python2
else
  python_command=python
fi

$python_command configure COPTFLAGS="-O3 -g" CXXOPTFLAGS="-O3 -g" FOPTFLAGS="-O3 -g" $ldflags_petsc  --with-cxx-dialect=C++11 $petsc_openmp --with-mpi-dir=$mpi_dir $configure_options --with-mumps-lib="$MUMPS_extra_lib" --prefix=$1/PETSC --with-debugging=0
make all test
make install

# if empty remove the folder
if [ ! "$(ls -A $1/PETSC)" ]; then
   rm -rf $1/PETSC
else
   #Mark the installation
   echo 1 > $1/PETSC/version
   exit 0
fi

