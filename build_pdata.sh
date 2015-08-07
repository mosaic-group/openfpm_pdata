#! /bin/bash

# Make a directory in /tmp/OpenFPM_data

echo "Directory: $1"
echo "Machine: $2"

mkdir /tmp/openfpm_pdata
mv * .[^.]* /tmp/openfpm_pdata
mv /tmp/openfpm_pdata OpenFPM_pdata

mkdir OpenFPM_pdata/src/config

git clone ssh://git@ppmcoremirror.dynu.com:2222/incardon/openfpm_vcluster.git OpenFPM_vcluster
git clone ssh://git@ppmcoremirror.dynu.com:2222/incardon/openfpm_devices.git OpenFPM_devices
git clone ssh://git@ppmcoremirror.dynu.com:2222/incardon/openfpm_data.git OpenFPM_data
git clone ssh://git@ppmcoremirror.dynu.com:2222/incardon/openfpm_io.git OpenFPM_IO
cd OpenFPM_data
git checkout develop
cd ..

cd "$1/OpenFPM_pdata"

if [ "$2" == "gin" ]
then
 echo "Compiling on gin\n"
 source ~/.bashrc
 sh ./autogen.sh
 module load gcc/4.9.2
 ./install
 make

 if [ $? -ne 0 ]; then exit 1 ; fi
 mpirun -np 2 ./src/pdata
 if [ $? -ne 0 ]; then exit 1 ; fi
 mpirun -np 3 ./src/pdata
 if [ $? -ne 0 ]; then exit 1 ; fi
 mpirun -np 4 ./src/pdata
 if [ $? -ne 0 ]; then exit 1 ; fi
 mpirun -np 5 ./src/pdata
 if [ $? -ne 0 ]; then exit 1 ; fi
 mpirun -np 6 ./src/pdata
 if [ $? -ne 0 ]; then exit 1 ; fi
 mpirun -np 7 ./src/pdata
 if [ $? -ne 0 ]; then exit 1 ; fi
 mpirun -np 8 ./src/pdata
 if [ $? -ne 0 ]; then exit 1 ; fi
 mpirun -np 9 ./src/pdata
 if [ $? -ne 0 ]; then exit 1 ; fi
 mpirun -np 10 ./src/pdata
 if [ $? -ne 0 ]; then exit 1 ; fi
 mpirun -np 11 ./src/pdata
 if [ $? -ne 0 ]; then exit 1 ; fi
 mpirun -np 12 ./src/pdata
 if [ $? -ne 0 ]; then exit 1 ; fi
elif [ "$2" == "wetcluster" ]
then
 echo "Compiling on wetcluster"

## produce the module path

 source ~/.bashrc
 module load gcc/4.9.2
 module load openmpi/1.8.1
 module load boost/1.54.0

 sh ./autogen.sh
 ./install --with-boost=/sw/apps/boost/1.54.0/  CXX=mpic++
 make
 if [ $? -ne 0 ]; then exit 1 ; fi

 ## Run on the cluster
 bsub -o output_run2.%J -K -n 2 -R "span[hosts=1]" "module load openmpi/1.8.1 ; module load gcc/4.9.2;  mpirun -np 2 ./src/pdata"
 if [ $? -ne 0 ]; then exit 1 ; fi
 bsub -o output_run3.%J -K -n 3 -R "span[hosts=1]" "module load openmpi/1.8.1 ; module load gcc/4.9.2;  mpirun -np 3 ./src/pdata"
 if [ $? -ne 0 ]; then exit 1 ; fi
 bsub -o output_run4.%J -K -n 4 -R "span[hosts=1]" "module load openmpi/1.8.1 ; module load gcc/4.9.2;  mpirun -np 4 ./src/pdata"
 if [ $? -ne 0 ]; then exit 1 ; fi
 bsub -o output_run5.%J -K -n 5 -R "span[hosts=1]" "module load openmpi/1.8.1 ; module load gcc/4.9.2;  mpirun -np 5 ./src/pdata"
 if [ $? -ne 0 ]; then exit 1 ; fi
 bsub -o output_run6.%J -K -n 6 -R "span[hosts=1]" "module load openmpi/1.8.1 ; module load gcc/4.9.2;  mpirun -np 6 ./src/pdata"
 if [ $? -ne 0 ]; then exit 1 ; fi
 bsub -o output_run7.%J -K -n 7 -R "span[hosts=1]" "module load openmpi/1.8.1 ; module load gcc/4.9.2;  mpirun -np 7 ./src/pdata"
 if [ $? -ne 0 ]; then exit 1 ; fi
 bsub -o output_run8.%J -K -n 8 -R "span[hosts=1]" "module load openmpi/1.8.1 ; module load gcc/4.9.2;  mpirun -np 8 ./src/pdata"
 if [ $? -ne 0 ]; then exit 1 ; fi
 bsub -o output_run9.%J -K -n 9 -R "span[hosts=1]" "module load openmpi/1.8.1 ; module load gcc/4.9.2;  mpirun -np 9 ./src/pdata"
 if [ $? -ne 0 ]; then exit 1 ; fi
 bsub -o output_run10.%J -K -n 10 -R "span[hosts=1]" "module load openmpi/1.8.1 ; module load gcc/4.9.2;  mpirun -np 10 ./src/pdata"
 if [ $? -ne 0 ]; then exit 1 ; fi
 bsub -o output_run11.%J -K -n 11 -R "span[hosts=1]" "module load openmpi/1.8.1 ; module load gcc/4.9.2;  mpirun -np 11 ./src/pdata"
 if [ $? -ne 0 ]; then exit 1 ; fi
 bsub -o output_run12.%J -K -n 12 -R "span[hosts=1]" "module load openmpi/1.8.1 ; module load gcc/4.9.2;  mpirun -np 12 ./src/pdata"
 if [ $? -ne 0 ]; then exit 1 ; fi
 # bsub -o output_run32.%J -K -n 32 "module load openmpi/1.8.1 ; module load gcc/4.9.2;  mpirun -np 32 ./src/vcluster"
 # if [ $? -ne 0 ]; then exit 1 ; fi
 # bsub -o output_run32.%J -K -n 64 "module load openmpi/1.8.1 ; module load gcc/4.9.2;  mpirun -np 64 ./src/vcluster"
 # if [ $? -ne 0 ]; then exit 1 ; fi
 # bsub -o output_run32.%J -K -n 128 "module load openmpi/1.8.1 ; module load gcc/4.9.2;  mpirun -np 128 ./src/vcluster"
 # if [ $? -ne 0 ]; then exit 1 ; fi

elif [ "$2" == "taurus" ]
then
 echo "Compiling on taurus"

 source /etc/profile
 echo "$PATH"
 module load gcc/4.8.2
 module load boost/1.55.0-gnu4.8
 module load openmpi/1.8.7
 module unload bullxmpi

 sh ./autogen.sh
 ./configure --with-metis=$METIS_ROOT CXX=mpic++
 make
 if [ $? -ne 0 ]; then exit 1 ; fi

# salloc --nodes=1 --ntasks-per-node=16 --time=00:05:00 --mem-per-cpu=1900 --partition=haswell bash -c "ulimit -s unlimited && mpirun -np 16 src/vcluster"
# if [ $? -ne 0 ]; then exit 1 ; fi
# salloc --nodes=2 --ntasks-per-node=16 --time=00:05:00 --mem-per-cpu=1900 --partition=haswell bash -c "ulimit -s unlimited && mpirun -np 32 src/vcluster"
# if [ $? -ne 0 ]; then exit 1 ; fi
# salloc --nodes=4 --ntasks-per-node=16 --time=00:05:00 --mem-per-cpu=1900 --partition=haswell bash -c "ulimit -s unlimited && mpirun -np 64 src/vcluster"
# if [ $? -ne 0 ]; then exit 1 ; fi
# salloc --nodes=8 --ntasks-per-node=16 --time=00:05:00 --mem-per-cpu=1900 --partition=haswell bash -c "ulimit -s unlimited && mpirun -np 128 src/vcluster"
# if [ $? -ne 0 ]; then exit 1 ; fi
# salloc --nodes=16 --ntasks-per-node=16 --time=00:5:00 --mem-per-cpu=1900 --partition=haswell bash -c "ulimit -s unlimited && mpirun -np 256 src/vcluster"
# if [ $? -ne 0 ]; then exit 1 ; fi

else
 echo "Compiling general"
 source ~/.bashrc
 sh ./autogen.sh
 ./install

 mpirun -np 2 ./src/pdata
 mpirun -np 3 ./src/pdata
 mpirun -np 4 ./src/pdata
fi



