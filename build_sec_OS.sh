#! /bin/bash

# Make a directory in /tmp/OpenFPM_pdata

echo "Directory: $1"
echo "Machine: $2"
echo "Branch: $BRANCH_NAME"

mkdir src/config


if [ "$2" == "windows10" ]; then
    echo "Compiling on windows10"

    echo "1" > input_install
    echo "2" >> input_install
    echo "y" >> input_install
    echo "1" >> input_install

    ./install -i "/home/jenkins/$BRANCH_NAME" < input_install

fi


