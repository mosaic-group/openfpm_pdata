#! /bin/bash

# Make a directory in /tmp/OpenFPM_pdata

echo "Directory: $1"
echo "Machine: $2"
echo "Branch: $6"

if [ "$2" == "gin" ]
then
 module load gcc/4.9.2
 if [ $? -ne 0 ]; then
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to complete the openfpm_pdata test \" }] }" https://hooks.slack.com/services/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
   exit 1 ;
 fi

 source $HOME/openfpm_vars_$6

 echo "$LD_LIBRARY_PATH"
 ldd ./src/pdata

 if [ x"$3" == x"no_test" ]; then
   exit 0;
 fi

 mpirun -np $3 ./src/pdata
 if [ $? -ne 0 ]; then 
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to complete the openfpm_pdata test \" }] }" https://hooks.slack.com/services/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
   exit 1 ;
 fi

elif [ "$2" == "taurus" ]
then

 source /etc/profile
 echo "$PATH"
 module load eigen/3.2.0
 module load suitesparse/4.2.1-gnu-multimkl
 module load boost/1.60.0
 module load gcc/5.3.0
 module load openmpi/1.10.2-gnu
 module unload bullxmpi
 
 export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/incard/PARMETIS/lib:/home/incard/METIS/lib:/home/incard/HDF5/lib"

 source $HOME/openfpm_vars_$6

 salloc --nodes=$4 --ntasks-per-node=$5 --time=00:15:00 --mem-per-cpu=1900 --partition=haswell bash -c "ulimit -s unlimited && mpirun -np $3 src/pdata --report_level=no"
 if [ $? -ne 0 ]; then
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to complete the openfpm_pdata test \" }] }" https://hooks.slack.com/services/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
   exit 1 ; 
 fi

else

 if [ x"$3" == x"no_test" ]; then
   exit 0;
 fi

 source $HOME/openfpm_vars_$6

 mpirun -np $3 ./src/pdata
 if [ $? -ne 0 ]; then
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to complete the openfpm_pdata test \" }] }" https://hooks.slack.com/services/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
   exit 1 ;
 fi
fi


