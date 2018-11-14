#! /bin/bash

# Make a directory in /tmp/OpenFPM_pdata

echo "Directory: $1"
echo "Machine: $2"
echo "Branch: $6"

exit 1

if [ "$2" == "gin" ]
then
 source "$HOME/.bashrc"
 module load gcc/4.9.2
 if [ $? -ne 0 ]; then
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$2 failed to complete the openfpm_pdata test \" }] }" https://hooks.slack.com/services/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
   exit 1 ;
 fi

 if [ x"$6" != x"" ]; then
   source $HOME/openfpm_vars_$6
 else
   source $HOME/openfpm_vars_master
 fi


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
 module load gcc/5.5.0
 module load openmpi/3.0.0-gnu5.5
 module unload bullxmpi

 if [ x"$6" != x"" ]; then
   source $HOME/openfpm_vars_$6
 else
   source $HOME/openfpm_vars_master
 fi

 salloc --nodes=$4 --ntasks-per-node=$5 --time=00:35:00 --mem-per-cpu=1900 --partition=haswell bash -c "ulimit -s unlimited && mpirun -np $3 src/pdata --report_level=no"
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


