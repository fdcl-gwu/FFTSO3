#ifndef _FDCL_FFTS2_MATRIX_HPP
#define _FDCL_FFTS2_MATRIX_HPP

#include <iostream>
#include <vector>
#include <math.h> // pow
#include <Eigen/Dense>

#ifndef _IMAGINARY_UNIT
#define _IMAGINARY_UNIT
const std::complex<double> I(0.0,1.0);    
#endif

using std::cout;
using std::endl;
using std::ostream;
using std::complex;

namespace fdcl
{
    template <class ScalarType> class FFTS2_matrix;
}

template <class ScalarType>
class fdcl::FFTS2_matrix
{
    public:
        int l_max;
        std::vector<Eigen::Matrix<ScalarType,Eigen::Dynamic,1> > M;
        FFTS2_matrix(){l_max=0;};
        ~FFTS2_matrix(){};
        FFTS2_matrix(int l_max); 

        void init(int l_max);
        Eigen::Matrix<ScalarType,Eigen::Dynamic,1>& operator[](int l); // return l-th matrix 
        ScalarType& operator()(int l, int m); // access (l,m)-th element of the l-th matrix
        void setRandom();
        void setZero();
        double norm();

        fdcl::FFTS2_matrix<double> real();  	

        template<typename _ScalarType>
            friend ostream& operator<<(ostream& os, const fdcl::FFTS2_matrix<_ScalarType>& M);  	

        fdcl::FFTS2_matrix<std::complex<double>> operator+(fdcl::FFTS2_matrix<std::complex<double>> const& M1);  	
        fdcl::FFTS2_matrix<ScalarType> operator+(fdcl::FFTS2_matrix<double> const& M2);  		
        fdcl::FFTS2_matrix<std::complex<double>> operator-(fdcl::FFTS2_matrix<std::complex<double>> const& M1);  	
        fdcl::FFTS2_matrix<ScalarType> operator-(fdcl::FFTS2_matrix<double> const& M2);  		

        fdcl::FFTS2_matrix<std::complex<double>> operator*(const std::complex<double>& c);  	
        fdcl::FFTS2_matrix<ScalarType> operator*(const double& c);  	

    private:
        void assert_index(int l);
        void assert_index(int l, int m);
};


namespace fdcl
{
    typedef fdcl::FFTS2_matrix<double> FFTS2_matrix_real;
    typedef fdcl::FFTS2_matrix<std::complex<double>> FFTS2_matrix_complex;
}

#endif
