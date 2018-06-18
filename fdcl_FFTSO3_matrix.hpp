// Fourier transform on SO(3) uses (2l+1)\times(2l+1) matrices for l\in{0,\ldots l_\max},
// whose row and column indicies are given by -l \leq m,n \leq l.
// This class provides a set of matrices with the dimensions described above. 
// The operator [l] is overloaded to yield the l-th matrix.
// The operator (l,m,n) is overloaded to access the element with the  -l \leq m,n \leq l

#ifndef _FDCL_FFTSO3_MATRIX_HPP
#define _FDCL_FFTSO3_MATRIX_HPP

#include <iostream>
#include <vector>
#include <math.h> // pow
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const complex<double> I(0.0,1.0);    

template <class ScalarType>
class fdcl_FFTSO3_matrix
{
public:
	int l_max;
	std::vector<Eigen::Matrix<ScalarType,Dynamic,Dynamic>> M;
	fdcl_FFTSO3_matrix(){};
	~fdcl_FFTSO3_matrix(){};
	fdcl_FFTSO3_matrix(int l_max); 
	void init(int l_max);
	Eigen::Matrix<ScalarType,Dynamic,Dynamic>& operator[](int l); // return l-th matrix 
	ScalarType& operator()(int l, int m, int n); // access (m,n)-th element of the l-th matrix
	void setRandom();
	void setZero();
	
	template<typename _ScalarType>
    friend ostream& operator<<(ostream& os, const fdcl_FFTSO3_matrix<_ScalarType>& M);  	

	fdcl_FFTSO3_matrix<complex<double>> operator+(fdcl_FFTSO3_matrix<complex<double>> const& M1);  	
	fdcl_FFTSO3_matrix<ScalarType> operator+(fdcl_FFTSO3_matrix<double> const& M2);  		

    fdcl_FFTSO3_matrix<complex<double>> operator*(const complex<double>& c);  	
    fdcl_FFTSO3_matrix<ScalarType> operator*(const double& c);  	

private:
	void assert_index(int l);
	void assert_index(int l, int m, int n);
};

typedef fdcl_FFTSO3_matrix<double> fdcl_FFTSO3_matrix_real;
typedef fdcl_FFTSO3_matrix<complex<double>> fdcl_FFTSO3_matrix_complex;

#endif