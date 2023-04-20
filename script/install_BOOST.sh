#!/bin/bash 

short_date=$(/bin/date +%m%d%y)
exec 2>>"install_BOOST$short_date.log"
set -x

source script/discover_os
discover_os

# check if the directory $1/BOOST exist

if [ -d "$1/BOOST" ]; then
  echo "BOOST already installed"
  exit 0
fi

rm boost_1_82_0.tar.bz2
wget http://ppmcore.mpi-cbg.de/upload/boost_1_82_0.tar.bz2
tar -xf boost_1_82_0.tar.bz2
cd boost_1_82_0
if [ x"$4" != x"" ]; then
	if [ -f $HOME/user-config.jam ]; then
		mv $HOME/user-config.jam $HOME/user-config.jam_bck
	fi
	if [ x"$5" != x"" ]; then
		echo "using gcc : $5.$6 : $4 ; " > $HOME/user-config.jam
	else
		echo "using gcc : : $4 ; " > $HOME/user-config.jam
	fi
fi

## When we are on powerPC we avoid booststrap for some reason will try xlcpp and most probably will fail we fall back to gcc

if [ x"$arch" == x"ppc64le" ]; then
	cd tools/build/src/engine
	./build.sh gcc
	cd ../../../../
	cp tools/build/src/engine/b2 .
	cp tools/build/src/engine/bjam .
else
	./bootstrap.sh --with-toolset=$3
fi


mkdir $1/BOOST
# Several flavours
if [ x"$platform" == x"osx" ]; then
    if [ x"$arch" == x"arm64" ]; then
        if [ x"$3" == x"" ]; then
            ./b2 -a -j $2 install --prefix=$1/BOOST address-model=64 architecture=arm abi=aapcs binary-format=mach-o toolset=clang  -sNO_LZMA=1 -sNO_ZSTD=1
        else
            ./b2 -a -j $2 install --prefix=$1/BOOST address-model=64 architecture=arm abi=aapcs binary-format=mach-o toolset=$3  -sNO_LZMA=1 -sNO_ZSTD=1
        fi
    else
        ./b2 -a -j $2 install --prefix=$1/BOOST address-model=64 architecture=x86 abi=sysv binary-format=mach-o toolset=clang  -sNO_LZMA=1 -sNO_ZSTD=1
    fi
else
    ./b2 -a -j $2 install --prefix=$1/BOOST  -sNO_LZMA=1 -sNO_ZSTD=1
fi

rm -rf boost_1_82_0

if [ -f $HOME/user-config.jam_bck ]; then
	mv $HOME/user-config.jam_bck $HOME/user-config.jam
fi
rm -rf boost_1_82_0.tar.bz2

