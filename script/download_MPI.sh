#!/bin/bash

# check if the directory $1/MPI exist

rm -rf openmpi-4.1.1
rm openmpi-4.1.1.tar.gz
wget http://ppmcore.mpi-cbg.de/upload/openmpi-4.1.1.tar.gz
tar -xvf openmpi-4.1.1.tar.gz

