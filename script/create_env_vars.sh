#! /bin/bash

prefix_dependencies="/usr/local"
prefix_openfpm="/usr/local"

if [[ "$OSTYPE" == "linux-gnu" -o "$OSTYPE" == "linux" -o "$OSTYPE" == "cygwin" ]]; then
  bash_library="export LD_LIBRARY_PATH=\"\$LD_LIBRARY_PATH:/$prefix_openfpm/openfpm_devices/lib:/$prefix_openfpm/openfpm_vcluster/lib"
else
  bash_library="export DYLD_LIBRARY_PATH=\"\$DYLD_LIBRARY_PATH:/$prefix_openfpm/openfpm_devices/lib:/$prefix_openfpm/openfpm_vcluster/lib"
fi

bash_path="export PATH=\""

if [ -d "$prefix_dependencies/MPI" ]; then
  bash_path="$bash_path$prefix_dependencies/MPI/bin:"
  bash_library="$bash_library:$prefix_dependencies/MPI/lib"
fi

if [ -d "$prefix_dependencies/METIS" ]; then
  bash_library="$bash_library:$prefix_dependencies/METIS/lib"
fi

if [ -d "$prefix_dependencies/PARMETIS" ]; then
  bash_library="$bash_library:$prefix_dependencies/PARMETIS/lib"
fi

if [ -d "$prefix_dependencies/BOOST" ]; then
  bash_library="$bash_library:$prefix_dependencies/BOOST/lib"
fi

if [ -d "$prefix_dependencies/HDF5" ]; then
  bash_library="$bash_library:$hdf5_lib"
fi

if [ -d "$prefix_dependencies/LIBHILBERT" ]; then
  bash_library="$bash_library:$prefix_dependencies/LIBHILBERT/lib"
fi

lin_alg_installed=""

if [ -d "$prefix_dependencies/PETSC" -a -f "$prefix_dependencies/PETSC/include/petsc.h" ]; then
  bash_library="$bash_library:$prefix_dependencies/PETSC/lib"
  lin_alg_installed="y"

if [ -d "$prefix_dependencies/OPENBLAS" ]; then
  bash_library="$bash_library:$prefix_dependencies/OPENBLAS/lib"
fi


if [ -d "$prefix_dependencies/SUITESPARSE"  -a -f "$prefix_dependencies/SUITESPARSE/include/umfpack.h" ]; then
  bash_library="$bash_library:$prefix_dependencies/SUITESPARSE/lib"
fi

# in cygwin we have to add to PATH additional directories
if [[ "$OSTYPE" == "cygwin" ]]; then
	bash_path="$bash_path:$prefix_dependencies/BOOST/bin:$prefix_dependencies/HDF5/bin"
fi

bash_path="$bash_path:\$PATH\""
bash_library="$bash_library\""

##### Writing openfpm_vars file

echo "$bash_path" > openfpm_vars
echo "$bash_library" >> openfpm_vars
echo "export PURE_PYTHON=1" >> openfpm_vars
