#! /bin/bash

# Make a directory in /tmp/OpenFPM_pdata

echo "Directory: $1"
echo "Machine: $2"

mkdir src/config

git submodule init
if [ $? -ne 0 ]; then
  echo -e "Configure\033[91;5;1m FAILED \033[0m"
  exit 1
fi

git submodule update
if [ $? -ne 0 ]; then
  echo -e "Configure\033[91;5;1m FAILED \033[0m"
  exit 1
fi

mkdir openfpm_numerics/src/config


if [ "$2" == "taurus" ]; then
    echo "Compiling on taurus with intel compiler"

    source /etc/profile
    echo "$PATH"
    module load gcc/4.9.3
    module load intel/2017.0.020

    echo "3" > input_install
    echo "1" >> input_install
    echo "1" >> input_install
    echo "24" >> input_install
    echo "y" >> input_install

    ./install -i "/scratch/p_ppm/openfpm_deps_intel" < input_install

fi


