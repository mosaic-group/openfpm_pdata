#!/bin/bash 

# check if the directory $1/BOOST exist
source script/discover_os
discover_os


if [ -d "$1/BOOST" ]; then
  echo "BOOST already installed"
  exit 0
fi

wget http://ppmcore.mpi-cbg.de/upload/boost_1_75_0.tar.bz2
if [ x"$platform" == x"msys" ]; then
#	bunzip2 boost_1_75_0.tar.bz2
#	tar -xvf boost_1_75_0.tar
#	cd boost_1_75_0
	pacman -Syu mingw-w64-x86_64-boost 
else
	tar -xvf boost_1_75_0.tar.bz2
	cd boost_1_75_0
	if [ x"$4" != x"" ]; then
		if [ -f "$HOME/user-config.jam" ]; then
			mv "$HOME/user-config.jam" "$HOME/user-config.jam_bck"
		fi
		if [ x"$5" != x"" ]; then
			echo "using gcc : $5.$6 : $4 ; " > "$HOME/user-config.jam"
		else
			echo "using gcc : : $4 ; " > "$HOME/user-config.jam"
		fi
	fi
	if [ -f "$HOME/user-config.jam_bck" ]; then
		mv "$HOME/user-config.jam_bck" "$HOME/user-config.jam"
	fi
fi
if [ x"$platform" == x"msys" ]; then
#	cd tools/build
#	./bootstrap.bat mingw
#	cd ..
#	cd ..
#	cp tools/build/b2.exe .
#	./b2 -j $2 install --prefix=$1/BOOST
	echo "BOOST installed"
else
	./bootstrap.sh --with-toolset=$3
	./b2 -j $2 install --prefix=$1/BOOST
fi
rm -rf boost_1_75_0

rm -rf boost_1_75_0.tar.bz2

