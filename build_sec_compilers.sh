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

source gcc/5.3.0


