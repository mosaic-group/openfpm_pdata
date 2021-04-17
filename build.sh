#! /bin/bash
# Make a directory in /tmp/OpenFPM_pdata

workspace=$1
hostname=$(hostname)
target=$3
comp_type=$4
branch=$5
with_gpu=$6

if [ x"$branch" == x"" ]; then
  echo "Getting branch from git"
  branch=$(git ls-remote --heads origin | grep $(git rev-parse HEAD) | cut -d / -f 3)
fi

echo "Directory: $workspace"
echo "Machine: $hostname"
echo "make target: $target"
echo "compilation type: $comp_type"
echo "Branch name: $branch"
echo "GPU compilation: $with_gpu"


if [ x"$hostname" == x"cifarm-centos-node.mpi-cbg.de"  ]; then
	echo "CentOS node"
	source /opt/rh/devtoolset-7/enable
fi

if [ x"$hostname" == x"cifarm-ubuntu-node"  ]; then
#	rm -rf $HOME/openfpm_dependencies/openfpm_pdata/$branch/
	echo "Ubuntu node"
	./install_MPI_mpich.sh $HOME/openfpm_dependencies/openfpm_pdata/$branch/ 4
	export PATH="/opt/bin:$PATH"
fi

if [ x"$hostname" == x"cifarm-mac-node.mpi-cbg.de"  ]; then
	echo "Mac node"
	export PATH="/usr/local/bin:$PATH"
#	rm -rf $HOME/openfpm_dependencies/openfpm_pdata/$branch/
fi

if [ x"$hostname" == x"falcon1" ]; then
#       rm -rf $HOME/openfpm_dependencies/openfpm_pdata/$branch/
        echo "falcon1 settings"
	if [ x"$comp_type" == x"intel" ]; then
        	module load parallel_studio_xe/2019u1
		dependency_dir=/projects/ppm/rundeck/openfpm_dependencies_intel/
        elif [ x"$with_gpu" == x"" ]; then
		mkdir /projects/ppm/rundeck/openfpm_dependencies_${branch}_no_cuda/
                dependency_dir=/projects/ppm/rundeck/openfpm_dependencies_${branch}_no_cuda/
	else
		mkdir /projects/ppm/rundeck/openfpm_dependencies_$branch/
		dependency_dir=/projects/ppm/rundeck/openfpm_dependencies_$branch/
	fi
else
	dependency_dir=$HOME/openfpm_dependencies/openfpm_pdata/$branch
	mkdir $HOME/openfpm_dependencies/openfpm_pdata/$branch
fi

if [ x"$with_gpu" == x"1" ]; then
	gpu_support=-g
else
	gpu_support=
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

if [ x"$comp_type" != x"full" ]; then
	installation_dir=" "
else
	installation_dir="--prefix=$HOME/openfpm_install/$branch"
fi

# force ssh to not use HostKey verification
#echo "StrictHostKeyChecking=no" > $HOME/.ssh/config
#chmod 600 $HOME/.ssh/config

foward_options="--enable-cuda-on-cpu"
install_options=
if [ x"$comp_type" == x"full" ]; then
        install_options="-s "
elif [ x"$comp_type" == x"intel" ]; then
        install_options="-s "
else
        install_options="-s -m "
fi

if [ x"$comp_type" == x"se_class" ]; then
	foward_options="$foward_options --enable-se-class1 --with-action-on-error=THROW_ON_ERROR"
elif [ x"$comp_type" == x"garbageinjv" ]; then
	foward_options="$foward_options  --enable-garbageinjv"
elif [ x"$comp_type" == x"asan" ]; then
        foward_options="$foward_options --enable-asan"
fi

echo "Installing with: ./install $gpu_support  -i $dependency_dir $install_options -c \"$installation_dir $foward_options  \"  "
./install $gpu_support -i $dependency_dir $install_options -c "$installation_dir $foward_options "
if [ $? -ne 0 ]; then
    echo "Fail to ./install"
    exit 1 ;
fi

# Check of we have to do a make install
if [ x"$comp_type" == x"full" ]; then
    mv $HOME/openfpm_vars $HOME/openfpm_vars_$branch
    make install
    if [ x"$?" != x"0" ]; then
        exit 1
    fi
else
    echo "Make install partial"
    if [ x"$comp_type" == x"intel" ]; then
        mv $HOME/openfpm_vars $HOME/openfpm_vars_intel
    else
        mv $HOME/openfpm_vars $HOME/openfpm_vars_$branch
    fi
    source $HOME/openfpm_vars_$branch
    if [ x"$hostname" == x"suitcase" ]; then
      echo "Running make on 1 cores"
      make VERBOSE=1 -j 1
    else
      make VERBOSE=1 -j 8
    fi
fi

if [ $? -ne 0 ]; then
   echo "Fail make install"
   exit 1 ;
fi

