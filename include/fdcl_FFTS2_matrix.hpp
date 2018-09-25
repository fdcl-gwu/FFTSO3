/** \file fdcl_FFTS2_matrix.hpp
 */
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

/** \brief Matrix class for spherical harmonics
 *
 * This class formulates a matrix class for harmonic analysis on the unit sphere.
 * Spherical harmonics is indexed by two integer variables, namely \f$(l,m)\f$ for \f$ 0\leq l\f$ and \f$ -l\leq m \leq l\f$. 
 * This class supports the following indexing:
 *   - \f$ F^l_{m} \in\mathbb{C} \f$ is accessed by \c F(l,m)
 *
 * Furthermore, several common operations, such as \c real(), \c setZero(), \c +, \c -, \c* are provided for convenience.
 *
 * This is defined as a template class, and they are implemented depending on the type of the variable as follows.
 *   - fdcl::FFTS2_matrix_complex
 *   - fdcl::FFTS2_matrix_real
 */
template <class ScalarType>
class fdcl::FFTS2_matrix
{
    public:
        int l_max;  /**< The maximum order for the index \f$ l \f$ */
        std::vector<Eigen::Matrix<ScalarType,Eigen::Dynamic,1> > M; /**< Buffer array */
        FFTS2_matrix(){l_max=0;};
        ~FFTS2_matrix(){};
        FFTS2_matrix(int l_max);   /**< Constructor with \c l_max */

        void init(int l_max); /**< Initialize buffer with \c l_max */
        Eigen::Matrix<ScalarType,Eigen::Dynamic,1>& operator[](int l);  /**< Operator to access the \f$ (2l+1) \f$ vector for the order \f$l\f$ */
        ScalarType& operator()(int l, int m); /**< Operator to access the \f$(l,m,n)\f$-th element */

        void setRandom(); /**< Randomize every element */
        void setZero(); /**< Set every element to zero */

        /** Find the sum of the matrix norm:
         * \f[
         * \sum_{l=0}^L \| F^l \|.
         * \f]
         */
        double norm();

        fdcl::FFTS2_matrix<double> real();  	 /**< Return the real parts */

        template<typename _ScalarType>
            friend ostream& operator<<(ostream& os, const fdcl::FFTS2_matrix<_ScalarType>& M);  	/**< Stream each matrix. For example, use it for \c std::out \c << \c F */  	

        fdcl::FFTS2_matrix<std::complex<double>> operator+(fdcl::FFTS2_matrix<std::complex<double>> const& M1);  /**< Operator overloading for addition */	
        fdcl::FFTS2_matrix<ScalarType> operator+(fdcl::FFTS2_matrix<double> const& M2);  	/**< Operator overloading for addition */	
        fdcl::FFTS2_matrix<std::complex<double>> operator-(fdcl::FFTS2_matrix<std::complex<double>> const& M1);  	/**< Operator overloading for subtraction */
        fdcl::FFTS2_matrix<ScalarType> operator-(fdcl::FFTS2_matrix<double> const& M2);  	/**< Operator overloading for subtraction */	

        fdcl::FFTS2_matrix<std::complex<double>> operator*(const std::complex<double>& c);  /**< Operator overloading for scalar multiplication */	
        fdcl::FFTS2_matrix<ScalarType> operator*(const double& c);  	/**< Operator overloading for scalar multiplication */	

    private:
        void assert_index(int l);
        void assert_index(int l, int m);
};


namespace fdcl
{
    typedef fdcl::FFTS2_matrix<double> FFTS2_matrix_real; /**< Template implementation of fdcl::FFTS2_matrix for real elements */ 
    typedef fdcl::FFTS2_matrix<std::complex<double>> FFTS2_matrix_complex; /**< Template implementation of fdcl::FFTS2_matrix for complex elements */ 
}

#endif
