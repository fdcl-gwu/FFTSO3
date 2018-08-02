#!/bin/bash

PYBIND_VER=2.2.2
INSTALL_DIR="/usr/local/include"
PYBIND_RELEASE_URL="https://github.com/pybind/pybind11/archive/v${PYBIND_VER}.tar.gz"
TEMP_DIR="$(mktemp -d)"

# This will download the latest eigen and install for the sy tem
echo "Going to install the latest Pybind11 to ${INSTALL_DIR}"
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

echo "We're going to download Pybind11 ${PYBIND_VER} and install to ${INSTALL_DIR}"
cd ${TEMP_DIR}
mkdir ${PYBIND_VER}
wget ${PYBIND_RELEASE_URL} -O ${TEMP_DIR}/${PYBIND_VER}.tar.gz
tar -xvzf ${PYBIND_VER}.tar.gz -C ./${PYBIND_VER} --strip-components=1

# move the eigen subdirectory to /usr/local/include
if [[ -d "${INSTALL_DIR}/pybind11" ]]; then
	echo "Pybind11 directory already  exists"
	read -p "Press Enter to remove ${INSTALL_DIR}/pybind11"
	rm -rf "${INSTALL_DIR}/pybind11"
fi

# copy to /usr/local/include
echo "Now copying to ${INSTALL_DIR}/pybind11"
sudo mv ${PYBIND_VER}/include/pybind11 ${INSTALL_DIR}

echo "Pybind11 ${PYBIND_VER} installed to ${INSTALL_DIR}/pybind11"
read -p "Press enter to exit"
