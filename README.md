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
* Complex/Real [Spherical Harmonics](https://en.wikipedia.org/wiki/Spherical_harmonics) on the Unit-Sphere
* Fast Fourier Transform on SO(3)
* Compute [Wigner D-matrix](https://en.wikipedia.org/wiki/Wigner_D-matrix)
* Compute [Clebsch-Gordon coefficients](https://en.wikipedia.org/wiki/Clebsch–Gordan_coefficients)


### Why should I use it?

* This packages utilizes [OpenMP](https://www.openmp.org) for accelerated computing for multicore processors.
* It implements Fast Fourier Transform (FFT) algorithms developed for SO(3).
* This is the only package that supports real harmonic analysis on SO(3).
* The computation in this package is verified by various unit-testing, utilizing [GoogleTest](https://github.com/google/googletest).
* For convenience, it supports the indexing consistent with mathematical expressions. More specifically, harmonics on SO(3) is indexed by three integers varying from negative values to positive ones:

<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=F^l_{m,n}\quad&space;\text{&space;for&space;}&space;0\leq&space;l<\infty,\;&space;-l\leq&space;m,n\leq&space;l." target="_blank"><img src="https://latex.codecogs.com/gif.latex?F^l_{m,n}\quad&space;\text{&space;for&space;}&space;0\leq&space;l<\infty,\;&space;-l\leq&space;m,n\leq&space;l." title="F^l_{m,n}\quad \text{ for } 0\leq l<\infty,\; -l\leq m,n\leq l." /></a>
</p>

* It is common that the indices `m,n` is mapped to non-negative values. 
In this package, the above element can be directly accessed by `F(l,m,n)` without any conversion. Also, the (2l+1) by (2l+1) matrix can be accessed by `F[l]`

* This package provides routines for [Clebsch-Gordon coefficients](https://en.wikipedia.org/wiki/Clebsch–Gordan_coefficients), or derivatives of harmonics that are not available elsewhere.
* It is based on the [Eigen library](http://eigen.tuxfamily.org/) supporting vectorization.

### What's the theoretical basis?

This library is based on

* T. Lee, "[Real Harmonic Analysis on the Special Orthogonal Group](https://arxiv.org/abs/1809.10533)," arXiv, 2018 (Real harmonic analysis on SO(3))

It also utilizes the results of 

* D. Varshalovich, A. Moskalaev, V. Khersonskii, "[Quantum Theory of Angular Momentum](https://www.amazon.com/Quantum-Theory-Angular-Momemtum-Varshalovich/dp/9971501074)," World Scientific, 1988 (Wigner-D function and spherical harmonics)
* G. Chirikjian, A. Kyatkin, "[Engineering Applications of Noncommutative Harmonic Analysis](https://www.amazon.com/Engineering-Applications-Noncommutative-Harmonic-Analysis/dp/0849307481)," CPC Press, 2000 (Operational properties)
* D. Marinucci, G. Peccati, "[Random Fields on the Sphere: Representation, Limit Theorems and Cosmological Applications](https://www.amazon.com/gp/product/0521175615/ref=oh_aui_search_detailpage?ie=UTF8&psc=1), Cambridge University Press, 2011 (Clebsch-Gordon coefficients)

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

**Notes for Windows:** While this package does not use any OS-specific command, it has not been tested in Windows. 

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

**Notes for Eigen library:** If the Eigen library is already installed, the command `git submodule update...` can be skipped. Instead, `CMakeList.txt` should be modified accordingly. See [Using Eigen in CMake Projets](https://eigen.tuxfamily.org/dox/TopicCMakeGuide.html).

## Examples

Four examples are provided. 

* `example0.cpp` minimal working example for real harmonic analysis on SO(3)
* `example1.cpp` elaborated example for real harmonic analysis on SO(3)
* `example2.cpp` elaborated example for complex harmonic analysis on SO(3)
* `example3.cpp` application to spherical image correlation introduced in this [article](https://arxiv.org/abs/1809.10533)


All of the example codes are under `FFTSO3/exmaple`, and once the above installation procedures are completed, the executable binary files are copied to the folder `FFTSO3/bin`

#### example0.cpp
This example illustrates how to perform fast forward transform of the trace function, and compute the inverse transform at the identity, using real harmonic analysis on SO(3).

```C++
#include <iostream>
#include "fdcl_FFTSO3.hpp"

// define a real-valued function on SO(3)
double func(Eigen::Matrix3d R)
{
    return R.trace();
}
    
int main()
{
    int l_max=2;  // the maximum order of Fourier transform
    fdcl::FFTSO3_real RFFTSO3(l_max); // FFTSO3_real object for real harmonic analysis on SO(3)
    fdcl::FFTSO3_matrix_real F(l_max); // FFTSO3_matrix_real object to save real-valued Fourier parameters
    Eigen::Matrix3d R; 
    double f=0.;

    F = RFFTSO3.forward_transform(func); // perform forward Fourier transform
    std::cout << "Fourier parameter F = " << std::endl << std::endl << F << std::endl; // show Fourier parameters

    R.setIdentity(); // R is set to the identity matrix
    f = RFFTSO3.inverse_transform(F,R); // compute the inverse transform at the identity
    cout << "f = " << f << std::endl; 

    return 0;
}
```
Excecute `FFTSO3/bin/example0` to get the followin results. 

```
Fourier parameter F =

l=0
-4.16334e-17

l=1
    0.333333 -4.53217e-17  7.28124e-17
-3.08183e-17     0.333333  8.25415e-18
-1.38104e-16 -5.13615e-18     0.333333

l=2
 5.48551e-18 -1.81924e-17  1.05714e-18  -3.3415e-17  3.14623e-18
 2.11998e-17 -4.16334e-17 -1.36215e-17  1.64077e-19  2.38253e-17
 -7.2492e-19 -9.99741e-18 -1.38778e-17  5.84227e-18  2.16787e-18
-3.84623e-17 -3.84837e-18  -1.0536e-17 -9.71445e-17  1.59126e-17
-8.54102e-19  4.87669e-17   7.8598e-18 -1.92921e-17  9.78959e-20


f = 3
```

For other examples, see the comments within the source file.

## User Manual

More detailed user manual is available at
[Doxygen Manual for Class Members (HTML)](https://fdcl-gwu.github.io/FFTSO3/doc/html/index.html)

## Relevant Projects
* [The SOFT Package:
FFTs on the Rotation Group](https://www.cs.dartmouth.edu/~geelong/soft/): complex harmonic analysis on SO(3)

## Contact
This library is developed by [Flight Dynamics and Control Lab](http://fdcl.seas.gwu.edu/) at The George Washington University, Washington DC. Contact [tylee@gwu.edu](mailto:tylee@gwu.edu) for question and comment.

## Acknowledgments 

This research has been supported in parts by NSF under the grant CMMI-1335008, and by AFOSR under the grant FA9550-18-1-0288.
