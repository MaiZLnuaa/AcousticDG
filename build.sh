#!/bin/sh

echo "******************************************";
echo "*     Discontinuous galerkin setup.      *";
echo "******************************************";
echo

echo "[0] Create some project directories if don't exist.";

mkdir -p 3rdParty
mkdir -p results

echo "[1] Get dependencies and external libraries.";

# Check if cmake is installed
which cmake > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "CMake found.";
else
    echo "CMake not found, installing...";
    sudo dnf install -y cmake
    echo "CMake installed.";
fi

# Check if g++ is installed
which g++ > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "G++ found.";
else
    echo "G++ not found, installing...";
    sudo dnf install -y gcc-c++
    export CC=gcc
    export CXX=g++
    echo "G++ installed.";
fi

# Check if gfortran is installed
which gfortran > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "Gfortran found.";
else
    echo "Gfortran not found, installing...";
    sudo dnf install -y gfortran
    echo "Gfortran installed.";
fi

# Check if Lapack/Blas is installed
dnf list installed blas-devel lapack-devel > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "Lapack/Blas found.";
else
    echo "Lapack/Blas not found, installing...";
    sudo dnf install blas-devel lapack-devel
    echo "Lapack/Blas installed.";
fi

# Check if libGLU is installed
dnf list installed mesa-libGLU > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "libGLU found.";
else
    echo "libGLU not found, installing...";
    sudo dnf install -y mesa-libGLU
    echo "libGLU installed.";
fi

# Check if libXft is installed
dnf list installed libXft > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "libXft found.";
else
    echo "libXft not found, installing...";
    sudo dnf install -y libXft
    echo "libXft installed.";
fi

# Check if libvtk9-dev is installed
dnf list installed vtk-devel > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "libvtk9-dev found.";
else
    echo "libvtk9-dev not found, installing...";
    sudo dnf install -y vtk-devel
    echo "libvtk9-dev installed.";
fi

# Check if FFTW3 is installed
dnf list installed fftw-devel > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "libfftw3-dev found.";
else
    echo "libfftw3-dev not found, installing...";
    sudo dnf install -y fftw3-devel
    echo "libfftw3-dev installed.";
fi

# Check if nlohmann-json is installed
dnf list installed json-devel > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "nlohmann-json3-dev found.";
else
    echo "nlohmann-json3-dev not found, installing...";
    sudo dnf install -y nlohmann-json-devel
    echo "nlohmann-json3-dev installed.";
fi

# No need to install GMSH or Eigen if already installed
echo "[2] GMSH and Eigen already installed, skipping installation.";

gmsh_version=4.13.1
# LINUX
if [ ! -d "3rdParty/gmsh" ]; then
	echo "Gmsh not found, installing...";
	# wget http://gmsh.info/bin/Linux/gmsh-${gmsh_version}-Linux64-sdk.tgz
	tar -xf gmsh-${gmsh_version}-Linux64-sdk.tgz
	rm -rf gmsh-${gmsh_version}-Linux64-sdk.tgz
	mv gmsh-${gmsh_version}-Linux64-sdk 3rdParty/gmsh
	echo "Gmsh installed."
else
	   echo "Gmsh found.";
fi
# macOS
# if [ ! -d "3rdParty/gmsh" ]; then
# 	echo "Gmsh not found, installing...";
# 	wget https://gmsh.info/bin/macOS/gmsh-${gmsh_version}-MacOSX-sdk.tgz
# 	tar -xf gmsh-${gmsh_version}-MacOSX-sdk.tgz
# 	rm -rf gmsh-${gmsh_version}-MacOSX-sdk.tgz
# 	mv gmsh-${gmsh_version}-MacOSX-sdk 3rdParty/gmsh
# 	echo "Gmsh installed."
# else
# 	echo "Gmsh found.";
# fi


cd 3rdParty/gmsh/
export FC=gfortran
export PATH=${PWD}/bin:${PWD}/lib:${PATH}
export INCLUDE=${PWD}/include:${INCLUDE}
export LIB=${PWD}/lib:${LIB}
export PYTHONPATH=${PWD}/lib:${PYTHONPATH} 
export DYLD_LIBRARY_PATH=${PWD}/lib:${DYLD_LIBRARY_PATH}
cd ../../

echo "[3] Build sources.";

rm -rf build/
mkdir build

cd build/
cmake ../ -DCMAKE_BUILD_TYPE=Release  -G "Unix Makefiles"
make -j
if [ $? -eq 0 ]; then
    echo "[end] Everything went successfully.";
else
    echo "[end] Error!";
fi
