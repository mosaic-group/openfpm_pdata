#! /bin/bash

# Make a directory in /tmp/OpenFPM_pdata

echo "Directory: $1"
echo "Machine: $2"
echo "Branch name: $5"

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



if [ "$2" == "gin" ]
then
 echo "Compiling on gin\n"

 source ~/.bashrc
 module load gcc/4.9.2
 mkdir $HOME/$5
 if [ x"$4" == x"full" ]; then
  ./install -i $HOME/$5  -s -c "--prefix=/home/jenkins/openfpm_install" < input_install
 elif [ x"$3" == x"numerics" ]; then
	 ./install -i $HOME/$5  -m -s -c "--prefix=/home/jenkins/openfpm_install" < input_install
  make $3
 else
  ./install -i $HOME/$5  -m -s -c "--prefix=/home/jenkins/openfpm_install --no-recursion" < input_install
  make $3
 fi
 if [ $? -ne 0 ]; then
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to complete the openfpm_pdata test \" }] }" https://hooks.slack.com/services/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
   exit 1 ;
 fi

 source $HOME/openfpm_vars

 if [ $? -ne 0 ]; then
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to complete the openfpm_pdata test \" }] }" https://hooks.slack.com/services/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
   exit 1 ; 
 fi

elif [ "$2" == "taurus" ]
then
 echo "Compiling on taurus"

 source /etc/profile
 echo "$PATH"
 module load eigen/3.2.0
 module load suitesparse/4.2.1-gnu-multimkl
 module load boost/1.60.0
 module load gcc/5.3.0
 module load openmpi/1.10.2-gnu
 module unload bullxmpi
 
 export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/incard/PARMETIS/lib:/home/incard/METIS/lib:/home/incard/HDF5/lib"

 mkdir /scratch/p_ppm/$5
 ./install -m -i "/scratch/p_ppm/$5" -s -c"CXX=mpic++ --no-recursion" < input_install
 make $3

 source $HOME/openfpm_vars

 if [ $? -ne 0 ]; then
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to complete the openfpm_pdata test \" }] }" https://hooks.slack.com/services/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
   exit 1 ; 
 fi
else
 echo "Compiling general"
 source ~/.bashrc

 mkdir $HOME/$5
 if [ x"$4" == x"full" ]; then
  ./install -i $HOME/$5  -s -c "--prefix=/Users/jenkins/openfpm_install" < input_install
 elif [ x"$3" == x"numerics" ]; then
  ./install -i $HOME/$5  -m -s -c "--prefix=/home/jenkins/openfpm_install" < input_install
  make $3
 else
  ./install -i $HOME/$5 -m -s -c "--prefix=/Users/jenkins/openfpm_install --no-recursion" < input_install
  make $3
 fi

 if [ $? -ne 0 ]; then
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to complete the openfpm_pdata test \" }] }" https://hooks.slack.com/services/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
   exit 1 ;
 fi

fi


