#! /bin/bash

# Make a directory in /tmp/OpenFPM_pdata

workspace=$1
hostname=$(hostname)
target=$3
comp_type=$4
branch=$5

if [ x"$branch" == x"" ]; then
  echo "Getting branch from git"
  branch=$(git ls-remote --heads origin | grep $(git rev-parse HEAD) | cut -d / -f 3)
fi

echo "Directory: $workspace"
echo "Machine: $hostname"
echo "make target: $target"
echo "compilation type: $comp_type"
echo "Branch name: $branch"


if [ x"$hostname" == x"cifarm-centos-node.mpi-cbg.de"  ]; then 
	./install_MPI_mpich.sh $HOME/openfpm_dependencies/openfpm_pdata/$branch/ 4
	echo 4 > $HOME/openfpm_dependencies/openfpm_pdata/$branch/MPI/version
fi

if [ x"$hostname" == x"cifarm-ubuntu-node.mpi-cbg.de"  ]; then
#	rm -rf $HOME/openfpm_dependencies/openfpm_pdata/$branch/
	echo "Continue"
fi

if [ x"$hostname" == x"cifarm-mac-node.mpi-cbg.de"  ]; then
	export PATH="/usr/local/bin:$PATH"
#	rm -rf $HOME/openfpm_dependencies/openfpm_pdata/$branch/
	echo "Continue"
fi

#### If you have a dep_dir file change the branch name to the dep_dir

dep_dir=$(cat dep_dir)
if [ x"$dep_dir" != x"" ]; then
  branch=$dep_dir
fi

mkdir src/config
mkdir openfpm_numerics/src/config

echo "Compiling general"

source ~/.bashrc
 
installation_dir=""
if [ x"$hostname" == x"sbalzarini-mac-15" ]; then
  installation_dir="--prefix=/Users/jenkins/openfpm_install"
else
  installation_dir="--prefix=$HOME/openfpm_install/$branch"
fi

# force ssh to not use HostKey verification
echo "StrictHostKeyChecking=no" > $HOME/.ssh/config
chmod 600 $HOME/.ssh/config

mkdir $HOME/openfpm_dependencies/openfpm_pdata/$branch
if [ x"$comp_type" == x"full" ]; then
  echo "Installing with: ./install -i $HOME/openfpm_dependencies/openfpm_pdata/$branch  -s -c \"$installation_dir\"  "
  ./install -i $HOME/openfpm_dependencies/openfpm_pdata/$branch  -s -c "$installation_dir"
  make install
  if [ $? -ne 0 ]; then
    curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$hostname failed to complete the openfpm_pdata test \" }] }" https://hooks.slack.com/services/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
    exit 1 ;
  fi
  mv $HOME/openfpm_vars $HOME/openfpm_vars_$branch
  source $HOME/openfpm_vars_$branch
elif [ x"$comp_type" == x"numerics" ]; then
  ./install -i $HOME/openfpm_dependencies/openfpm_pdata/$branch  -m -s -c "$installation_dir"

  if [ $? -ne 0 ]; then
    curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$hostname failed to complete the openfpm_pdata test \" }] }" https://hooks.slack.com/services/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
    exit 1 ;
  fi
  mv $HOME/openfpm_vars $HOME/openfpm_vars_$branch
  source $HOME/openfpm_vars_$branch
  make VERBOSE=1  -j 8
else
  echo "Installing with: ./install -i $HOME/openfpm_dependencies/openfpm_pdata/$branch -m -s -c \"$installation_dir --no-recursion\""
  ./install -i $HOME/openfpm_dependencies/openfpm_pdata/$branch -m -s -c "$installation_dir"

  if [ $? -ne 0 ]; then
    curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$hostname failed to complete the openfpm_pdata test \" }] }" https://hooks.slack.com/services/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
    exit 1 ;
  fi
  mv $HOME/openfpm_vars $HOME/openfpm_vars_$branch
  source $HOME/openfpm_vars_$branch

  echo "------------- DEBUGGING ------------------"
  echo "$PATH"
  echo "-------------------------------------------"
  make VERBOSE=1 -j 8
fi

if [ $? -ne 0 ]; then
   curl -X POST --data "payload={\"icon_emoji\": \":jenkins:\", \"username\": \"jenkins\"  , \"attachments\":[{ \"title\":\"Error:\", \"color\": \"#FF0000\", \"text\":\"$hostname failed to complete the openfpm_pdata test \" }] }" https://hooks.slack.com/services/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
   exit 1 ;
fi



