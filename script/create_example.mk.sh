#! /bin/bash

prefix_dependencies="/usr/local"
prefix_openfpm="/usr/local"

openmp_libs="$(cat openmp_libs)"
openmp_flags="$(cat openmp_flags)"
cuda_include_dirs="$(cat cuda_include)"
cuda_on_cpu="$(cat cuda_on_cpu)"
cuda_lib="$(cat cuda_lib)"
cuda_options="$(cat cuda_options)"

if [ -d "$prefix_dependencies/HDF5/lib" ]; then
  hdf5_lib=$prefix_dependencies/HDF5/lib
  hdf5_lib_dir=-L$prefix_dependencies/HDF5/lib
elif [ -d "$prefix_dependencies/HDF5/lib64" ]; then
  hdf5_lib=$prefix_dependencies/HDF5/lib64
  hdf5_lib_dir=-L$prefix_dependencies/HDF5/lib64
fi

lin_alg_dir=""
lin_alg_lib=""
lin_alg_inc=""

if [ -d "$prefix_dependencies/PETSC" ]; then
    lin_alg_dir="$lin_alg_dir -L$prefix_dependencies/PETSC/lib"
    lin_alg_lib="$lin_alg_lib -lpetsc"
    lin_alg_inc="$lin_alg_inc -I$prefix_dependencies/PETSC/include"
fi

if [ -d "$prefix_dependencies/OPENBLAS" ]; then
    lin_alg_dir="$lin_alg_dir -L$prefix_dependencies/OPENBLAS/lib"
    lin_alg_lib="$lin_alg_lib -lopenblas"
    lin_alg_inc="$lin_alg_inc -I$prefix_dependencies/OPENBLAS/include"
fi

if [ -d "$prefix_dependencies/SUITESPARSE"  -a -f "$prefix_dependencies/SUITESPARSE/include/umfpack.h" ]; then
    lin_alg_dir="$lin_alg_dir -L$prefix_dependencies/SUITESPARSE/lib"
    lin_alg_inc="$lin_alg_inc -I$prefix_dependencies/SUITESPARSE/include"
    lin_alg_lib="$lin_alg_lib -lumfpack -lamd -lbtf -lcamd -lccolamd -lcholmod -lcolamd -lcxsparse -lklu -ldl -lrbio -lspqr -lsuitesparseconfig"
fi

echo "INCLUDE_PATH= $cuda_include_dirs $openmp_flags  -I.  -I$prefix_openfpm/openfpm_numerics/include -I$prefix_openfpm/openfpm_pdata/include/config -I$prefix_openfpm/openfpm_pdata/include -I$prefix_openfpm/openfpm_data/include -I$prefix_openfpm/openfpm_vcluster/include -I$prefix_openfpm/openfpm_io/include -I$prefix_openfpm/openfpm_devices/include -I$prefix_dependencies/VCDEVEL/include  -I$prefix_dependencies/METIS/include -I$prefix_dependencies/PARMETIS/include -I$prefix_dependencies/BOOST/include -I$prefix_dependencies/HDF5/include -I$prefix_dependencies/LIBHILBERT/include  $lin_alg_inc -I$prefix_dependencies/BLITZ/include -I$prefix_dependencies/ALGOIM/include  -I$prefix_dependencies/SUITESPARSE/include " > example.mk
echo "LIBS_PATH=-L$prefix_openfpm/openfpm_devices/lib -L$prefix_openfpm/openfpm_pdata/lib  -L$prefix_openfpm/openfpm_vcluster/lib -L$prefix_dependencies/VCDEVEL/lib  -L$prefix_dependencies/METIS/lib -L$prefix_dependencies/PARMETIS/lib  -L$prefix_dependencies/BOOST/lib $hdf5_lib_dir -L$prefix_dependencies/LIBHILBERT/lib  $lin_alg_dir " >> example.mk
if [ x"$cuda_on_cpu" == x"YES" ]; then
   echo "CUDA_ON_CPU=YES" >> example.mk
fi
echo "LIBS=$openmp_flags $openmp_libs -lvcluster -lofpm_pdata -lofpmmemory -lparmetis -lmetis -lboost_iostreams -lboost_program_options -lhdf5 -llibhilbert -lVc  $cuda_lib $lin_alg_lib -ldl -lboost_filesystem -lboost_system $optional_boost" >> example.mk
echo "LIBS_NVCC=-Xcompiler=$openmp_flags -lvcluster -lofpm_pdata -lofpmmemory -lparmetis -lmetis -lboost_iostreams -lboost_program_options -lhdf5 -llibhilbert -lVc  $cuda_lib $lin_alg_lib -ldl -lboost_filesystem -lboost_system $optional_boost" >> example.mk
echo "INCLUDE_PATH_NVCC=-Xcompiler=$openmp_flags "$cuda_options"  -I. -I$prefix_openfpm/openfpm_numerics/include -I$prefix_openfpm/openfpm_pdata/include/config -I$prefix_openfpm/openfpm_pdata/include -I$prefix_openfpm/openfpm_data/include -I$prefix_openfpm/openfpm_vcluster/include -I$prefix_openfpm/openfpm_io/include -I$prefix_openfpm/openfpm_devices/include -I$prefix_dependencies/METIS/include -I$prefix_dependencies/PARMETIS/include -I$prefix_dependencies/BOOST/include -I$prefix_dependencies/HDF5/include -I$prefix_dependencies/LIBHILBERT/include  $lin_alg_inc -I$prefix_dependencies/BLITZ/include -I$prefix_dependencies/ALGOIM/include  -I$prefix_dependencies/SUITESPARSE/include " >> example.mk
