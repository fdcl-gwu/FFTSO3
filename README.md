[![Build Status](https://travis-ci.org/fdcl-gwu/FFTSO3.svg?branch=master)](https://travis-ci.org/fdcl-gwu/FFTSO3)

# FAST FOURIER TRANSFORM ON SO(3)

This package provides various tools for harmonic analysis on the special orthogonal group, SO(3), and also for spherical harmonics on the unit sphere. 

## Introduction 

### What does it do?

* Complex Noncommutative Harmonic Analysis on SO(3), (i.e., fast Forward transform and inverse transform)
* Real Noncommutative Harmonic Analysis on SO(3)
* Complex Spherical Harmonics
* Real Spherical Harmonics 


### Why should I use it?

* This packages utilizes [OpenMP](https://www.openmp.org) for accelerated computing with multithreaded computing
* It implements Fast Fourier Transform (FFT) algorithms developed for SO(3)
* This is the only package that supports real harmonic analysis on SO(3)
* The computation in this package is verified by various unit-testing, utilizing [GoogleTeset](https://github.com/google/googletest)
* For convenience, it supports the indexing consistent with mathematical expressions. More specifically, harmonics on SO(3) is indexed by three integers varying from negative values to positive ones:

	<center><a href="https://www.codecogs.com/eqnedit.php?latex=F^{l}_{m,n}&space;\text{&space;for&space;}&space;0\leq&space;l&space;\text{&space;and&space;}&space;-l\leq&space;m,n\leq&space;l" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F^{l}_{m,n}&space;\text{&space;for&space;}&space;0\leq&space;l&space;\text{&space;and&space;}&space;-l\leq&space;m,n\leq&space;l" title="F^{l}_{m,n} \text{ for } 0\leq l \text{ and } -l\leq m,n\leq l" /></a></center>
		
	It is common that the indicies `m,n` is mapped to non-negative values. 
	In this package, the above element can be direclty accessed by `F(l,m,n)` without any conversion. 
	Also, the (2l+1) by (2l+1) matrix can be accessed by `F[l]`


## Installation

### Required Libraries
The following libraries are required:

* [git](https://git-scm.com)
* [CMake](https://cmake.org)
* [OpenMP](https://www.openmp.org) 
* [Eigen](http://eigen.tuxfamily.org/)

Please follow the instruction at each link to install the libraries. 
The Eigen library can be skipped, as it is included as a git submodule in this package.  

**Tips for Mac:**  macOS comes with git. For other libraries, it is first recommended to install [Homebrew](https://brew.sh), a software package manager for Mac. Then, the above libraries can be installed by

	brew install cmake
	brew install eigen
	brew install llvm

**Notes for Windows:** While this pacakge does not use any OS-specific command, it has not been tested in Windows

### Installation 

Execute the following commands at a folder above this package will be installed:

	git clone https://github.com/fdcl-gwu/FFTSO3.git
	cd FFTSO3/build
	cmake ..
	make
	../bin/fftso3_unit_test

The last command executes unit-testing, the installation is succesful if it prints out the following message at the end

	...
	...
	[ PASSED ] 12 tests.



## Example

## User Manual

[Doxygen Manual for Class Members](https://fdcl-gwu.github.io/FFTSO3/doc/html/index.html)

