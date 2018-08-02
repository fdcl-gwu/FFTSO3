#!/bin/bash

EIGEN_VER=3.3.4
INSTALL_DIR="/usr/local/include"
EIGEN_RELEASE_URL="https://github.com/eigenteam/eigen-git-mirror/archive/${EIGEN_VER}.tar.gz"
TEMP_DIR="$(mktemp -d)"

# This will download the latest eigen and install for the sy tem
echo "Going to install the latest Eigen"
read -p "Press Enter to continue........"

# download  latest release of eigen
if [[ ! "$TEMP_DIR" || ! -d "$TEMP_DIR" ]]; then
	echo "Could not create temp dir"
	exit 1
fi

# delete the temp directory on cleanup
function cleanup {
    rm -rf "$TEMP_DIR"
    echo "Deleted temp working directory $TEMP_DIR"
}

trap cleanup EXIT

echo "We're going to download Eigen ${EIGEN_VER} and install to ${INSTALL_DIR}"
cd ${TEMP_DIR}
mkdir ${EIGEN_VER}
wget ${EIGEN_RELEASE_URL} -O ${TEMP_DIR}/${EIGEN_VER}.tar.gz
tar -xvzf ${EIGEN_VER}.tar.gz -C ./${EIGEN_VER} --strip-components=1

# move the eigen subdirectory to /usr/local/include
if [[ -d "${INSTALL_DIR}/Eigen" ]]; then
	echo "Eigen directory already  exists"
	read -p "Press Enter to remove ${INSTALL_DIR}/Eigen"
	rm -rf "${INSTALL_DIR}/Eigen"
fi

# delete the unsupported directory
if [[ -d "${INSTALL_DIR}/unsupported" ]]; then
    echo "Eigen unsupported directory already exists"
    read -p "Press Enter to remove ${INSTALL_DIR}/unsupported"
    rf -rf "${INSTALL_DIR}/Eigen"
fi

# echo "Now going to copy Eigen headers to /usr/include/local"
# read -p "Press enter to continue"

# # copy to /usr/local/include
# echo "Now copying to ${INSTALL_DIR}/Eigen"
# sudo mv ${EIGEN_VER}/Eigen ${INSTALL_DIR}

# # copy to usr/local/include
# echo "Now copying to ${INSTALL_DIR}/unsupported"
# sudo mv ${EIGEN_VER}/unsupported ${INSTALL_DIR}

# echo "Eigen ${EIGEN_VER} installed to ${INSTALL_DIR}/Eigen and ${INSTALL_DIR}/unsupported"

echo "Going to install Eigen using CMake"
read -p "Press enter to continue"
cd ${EIGEN_VER}
mkdir build
cd build
cmake ..
sudo make install

read -p "Press enter to exit"
