#! /bin/bash

gpu_support=1
prefix_dependencies="/usr/local"
prefix_openfpm="/usr/local"
configure_options=""

configure_options="CXX=mpic++ $configure_options"

if [ -d "$prefix_dependencies/PETSC" ]; then
    configure_options="$configure_options --with-petsc=$prefix_dependencies/PETSC "
fi

if [ -d "$prefix_dependencies/BOOST" ]; then
    configure_options=" $configure_options --with-boost=$prefix_dependencies/BOOST "
fi
if [ -d "$prefix_dependencies/HDF5" ]; then
    configure_options=" $configure_options --with-hdf5=$prefix_dependencies/HDF5/ "
fi
if [ -d "$prefix_dependencies/LIBHILBERT" ]; then
    configure_options=" $configure_options --with-libhilbert=$prefix_dependencies/LIBHILBERT "
fi

if [ -d "$prefix_dependencies/BLITZ" ]; then
    configure_options=" $configure_options --with-blitz=$prefix_dependencies/BLITZ "
fi

if [ -d "$prefix_dependencies/ALGOIM" ]; then
    configure_options=" $configure_options --with-algoim=$prefix_dependencies/ALGOIM "
fi

if [ -d "$prefix_dependencies/PARMETIS" ]; then
    configure_options=" $configure_options --with-parmetis=$prefix_dependencies/PARMETIS "
fi

if [ -d "$prefix_dependencies/METIS" ]; then
	configure_options=" $configure_options --with-metis=$prefix_dependencies/METIS "
fi

if [ -d "$prefix_dependencies/VCDEVEL" ]; then
    configure_options=" $configure_options --with-vcdevel=$prefix_dependencies/VCDEVEL "
fi

if [ -d "$prefix_dependencies/OPENBLAS" ]; then
    configure_options=" $configure_options --with-blas=$prefix_dependencies/OPENBLAS/"
fi

if [ -d "$prefix_dependencies/SUITESPARSE"  -a -f "$prefix_dependencies/SUITESPARSE/include/umfpack.h" ]; then
    configure_options="$configure_options --with-suitesparse=$prefix_dependencies/SUITESPARSE "
fi

if [ -d "$prefix_dependencies/EIGEN" ]; then
    configure_options=" $configure_options --with-eigen=$prefix_dependencies/EIGEN "
    lin_alg_inc="$lin_alg_inc -I$prefix_dependencies/EIGEN"
fi

if [[ "$OSTYPE" == "linux-gnu" -o "$OSTYPE" == "linux" ]]; then
    lin_alg_lib="$lin_alg_lib -lrt"
fi

if [ x"$gpu_support" == x"1" ]; then
  configure_options=" $configure_options --with-cuda-on-backend=CUDA"
fi

./configure $$configure_options

