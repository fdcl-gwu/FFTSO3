#ifndef _FDCL_FFTS2_MATRIX_HPP
#define _FDCL_FFTS2_MATRIX_HPP

#include <iostream>
#include <vector>
#include <math.h> // pow
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#ifndef _IMAGINARY_UNIT
#define _IMAGINARY_UNIT
const complex<double> I(0.0,1.0);    
#endif

template <class ScalarType>
class fdcl_FFTS2_matrix
{
public:
	int l_max;
	std::vector<Eigen::Matrix<ScalarType,Dynamic,1>> M;
	fdcl_FFTS2_matrix(){l_max=0;};
	~fdcl_FFTS2_matrix(){};
	fdcl_FFTS2_matrix(int l_max); 

	void init(int l_max);
	Eigen::Matrix<ScalarType,Dynamic,1>& operator[](int l); // return l-th matrix 
	ScalarType& operator()(int l, int m); // access (l,m)-th element of the l-th matrix
	void setRandom();
	void setZero();
    double norm();
	
	fdcl_FFTS2_matrix<double> real();  	

	template<typename _ScalarType>
    friend ostream& operator<<(ostream& os, const fdcl_FFTS2_matrix<_ScalarType>& M);  	

	fdcl_FFTS2_matrix<complex<double>> operator+(fdcl_FFTS2_matrix<complex<double>> const& M1);  	
	fdcl_FFTS2_matrix<ScalarType> operator+(fdcl_FFTS2_matrix<double> const& M2);  		
	fdcl_FFTS2_matrix<complex<double>> operator-(fdcl_FFTS2_matrix<complex<double>> const& M1);  	
	fdcl_FFTS2_matrix<ScalarType> operator-(fdcl_FFTS2_matrix<double> const& M2);  		
    
    fdcl_FFTS2_matrix<complex<double>> operator*(const complex<double>& c);  	
    fdcl_FFTS2_matrix<ScalarType> operator*(const double& c);  	

private:
	void assert_index(int l);
	void assert_index(int l, int m);
};

typedef fdcl_FFTS2_matrix<double> fdcl_FFTS2_matrix_real;
typedef fdcl_FFTS2_matrix<complex<double>> fdcl_FFTS2_matrix_complex;

#endif
