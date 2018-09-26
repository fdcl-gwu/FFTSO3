[![Build Status](https://travis-ci.org/fdcl-gwu/FFTSO3.svg?branch=master)](https://travis-ci.org/fdcl-gwu/FFTSO3)

# FAST FOURIER TRANSFORM ON SO(3)

This C++ software package provides various tools for harmonic analysis on the special orthogonal group, 

<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\mathrm{SO(3)}=\{&space;R\in\Re^{3\times&space;3}\,|\,&space;R^TR=I_{3\times&space;3},\quad&space;\mathrm{det}[R]=&plus;1\}," target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathrm{SO(3)}=\{&space;R\in\Re^{3\times&space;3}\,|\,&space;R^TR=I_{3\times&space;3},\quad&space;\mathrm{det}[R]=&plus;1\}," title="\mathrm{SO(3)}=\{ R\in\Re^{3\times 3}\,|\, R^TR=I_{3\times 3},\quad \mathrm{det}[R]=+1\}," /></a>
</p>

which is the configuration space for the attitude dynamics of a rigid body.


## Introduction 

### What does it do?

* Complex/Real [Noncommutative Harmonic Analysis](https://en.wikipedia.org/wiki/Noncommutative_harmonic_analysis) on SO(3)
* Complex/Real [Spherical Harmonics](https://en.wikipedia.org/wiki/Spherical_harmonics)
* Fast Fourier Transform


### Why should I use it?

* This packages utilizes [OpenMP](https://www.openmp.org) for accelerated computing for multicore processors.
* It implements Fast Fourier Transform (FFT) algorithms developed for SO(3).
* This is the only package that supports real harmonic analysis on SO(3).
* The computation in this package is verified by various unit-testing, utilizing [GoogleTeset](https://github.com/google/googletest).
* For convenience, it supports the indexing consistent with mathematical expressions. More specifically, harmonics on SO(3) is indexed by three integers varying from negative values to positive ones:

<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=F^l_{m,n}\quad&space;\text{&space;for&space;}&space;0\leq&space;l<\infty,\;&space;-l\leq&space;m,n\leq&space;l." target="_blank"><img src="https://latex.codecogs.com/gif.latex?F^l_{m,n}\quad&space;\text{&space;for&space;}&space;0\leq&space;l<\infty,\;&space;-l\leq&space;m,n\leq&space;l." title="F^l_{m,n}\quad \text{ for } 0\leq l<\infty,\; -l\leq m,n\leq l." /></a>
</p>

	It is common that the indicies `m,n` is mapped to non-negative values. 
	In this package, the above element can be direclty accessed by `F(l,m,n)` without any conversion. Also, the (2l+1) by (2l+1) matrix can be accessed by `F[l]`

* This package provides routines for [Clebsch-Gordon coefficients](https://en.wikipedia.org/wiki/Clebsch–Gordan_coefficients), or derivatives of harmonics that are not available elsewhere.
* It is based on the [Eigen library](http://eigen.tuxfamily.org/) supporting vectorization.

## Installation

### Required Libraries
The following libraries are required:

* [git](https://git-scm.com)
* [CMake](https://cmake.org)
* [OpenMP](https://www.openmp.org) 
* [(Eigen)](http://eigen.tuxfamily.org/)

Please follow the instruction at each link to install the libraries, except the Eigen library that is included as a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) in this repository.

**Tips for Mac:**  First, install [XCode Command Line Tools](http://railsapps.github.io/xcode-command-line-tools.html) for gcc compilers and git. For other libraries, it is recommended to install [Homebrew](https://brew.sh), a software package manager for Mac. Then, the CMake and OpenMP can be installed by

```
	brew install cmake
	brew install llvm	
```
	
**Tips for Linux:** For Debian linux, such as Ubuntu, the above packages can be installed by
 
```
	sudo apt-get install git-core 
	sudo apt-get install cmake
	sudo apt-get install libomp-dev
```

**Notes for Windows:** While this pacakge does not use any OS-specific command, it has not been tested in Windows. 

### Compile 

Execute the following commands at a folder above this package will be installed:

```
	git clone https://github.com/fdcl-gwu/FFTSO3.git
	cd FFTSO3
	git submodule update --init --recursive
	cd build
	cmake ..
	make
	../bin/fftso3_unit_test
```

The last command executes unit-testing, and the installation is succesful if it prints out the following message at the end

```
	...
	...
	[ PASSED ] 12 tests.
```

**Notes for Eigen library:** If the Eigen library is alraedy installed, the command `git submodule update...` can be skipped. Instead, `CMakeList.txt` should be modified accordingly. See [Using Eigen in CMake Projets](https://eigen.tuxfamily.org/dox/TopicCMakeGuide.html).

## Example

## User Manual

[Doxygen Manual for Class Members](https://fdcl-gwu.github.io/FFTSO3/doc/html/index.html)

