#!/bin/bash

# 3.10.2 has an issue with the Finding boost 1.66
VERSION=3.11
BUILD="3"
TEMP_DIR="$(mktemp -d)"

echo "This will download and install the latest stable version of CMAKE"
echo "First we'll remove cmake"

read -p "Press Enter to continue"

sudo apt-get -y purge cmake
sudo apt-get -y install checkinstall

echo "Now going to download cmake v$VERSION.$BUILD"

read -p "Press Enter to continue"

if [[ ! "$TEMP_DIR" || ! -d "$TEMP_DIR" ]]; then
	echo "Could not create temp dir"
	exit 1
fi

cd ${TEMP_DIR}
wget https://cmake.org/files/v$VERSION/cmake-$VERSION.$BUILD.tar.gz
tar -xzvf cmake-$VERSION.$BUILD.tar.gz
cd cmake-$VERSION.$BUILD/

echo "Now going to configure cmake"
read -p "Press Enter to continue"

./bootstrap

echo "Now build cmake"
make -j 4
sudo checkinstall make install
