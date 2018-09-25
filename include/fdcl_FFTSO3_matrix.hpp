/** \file fdcl_FFTSO3_matrix.hpp
 */

#ifndef _FDCL_FFTSO3_MATRIX_HPP
#define _FDCL_FFTSO3_MATRIX_HPP

#include <iostream>
#include <vector>
#include <math.h> // pow
#include <Eigen/Dense>

using std::cout;
using std::endl;
using std::ostream;
using std::complex;
using std::max;
using std::min;

#ifndef _IMAGINARY_UNIT
#define _IMAGINARY_UNIT
const complex<double> I(0.0,1.0);    
#endif

namespace fdcl
{
    template <class ScalarType> class FFTSO3_matrix;
}

/** \brief Matrix class for harmonic analysis on SO(3)
 *
 * This class formulates a matrix class for harmonic analysis on SO(3). 
 * Harmonics on SO(3) is indexed by three integer variables, namely \f$(l,m,n)\f$ for \f$ 0\leq l\f$ and \f$ -l\leq m,n \leq l\f$. 
 * This class supports the following indexing:
 *   - \f$ F^l_{m,n} \in\mathbb{C} \f$ is accessed by \c F(l,m,n)
 *   - \f$ F^l \in \mathbb{C}^{(2l+1)\times(2l+1)}\f$ is accessed by \c F[l]
 *
 * Furthermore, several common operations, such as \c real(), \c setZero(), \c +, \c -, \c* are provided for convenience.
 *
 * This is defined as a template class, and they are implemented depending on the type of the variable as follows.
 *   - fdcl::FFTSO3_matrix_complex
 *   - fdcl::FFTSO3_matrix_real
 */
template <class ScalarType>
class fdcl::FFTSO3_matrix
{
public:
	int l_max; /**< The maximum order for the index \f$ l \f$ */
	std::vector<Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic>> M; /**< Buffer array */

	FFTSO3_matrix(){};
	~FFTSO3_matrix(){};
	FFTSO3_matrix(int l_max);  /**< Constructor with \c l_max */

	void init(int l_max); /**< Initialize buffer with \c l_max */

	Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic>& operator[](int l); /**< Operator to access the \f$ (2l+1)\times(2l+1) \f$ matrix for the order \f$l\f$ */
	ScalarType& operator()(int l, int m, int n); /**< Operator to access the \f$(l,m,n)\f$-th element */

	void setRandom(); /**< Randomize every element */
	void setZero(); /**< Set every element to zero */

    /** Find the sum of the matrix norm:
     * \f[
     * \sum_{l=0}^L \| F^l \|.
     * \f]
     */
    double norm(); 
	
	fdcl::FFTSO3_matrix<double> real(); /**< Return the real parts */

	template<typename _ScalarType>
    friend ostream& operator<<(ostream& os, const fdcl::FFTSO3_matrix<_ScalarType>& M);  	/**< Stream each matrix. For example, use it for \c std::out \c << \c F */

	fdcl::FFTSO3_matrix<complex<double>> operator+(fdcl::FFTSO3_matrix<complex<double>> const& M1); /**< Operator overloading for addition of two fdcl::FFTSO3_matrix objects */
	fdcl::FFTSO3_matrix<ScalarType> operator+(fdcl::FFTSO3_matrix<double> const& M2); /**< Operator overloading for addition of two fdcl::FFTSO3_matrix objects with real types */
	fdcl::FFTSO3_matrix<complex<double>> operator-(fdcl::FFTSO3_matrix<complex<double>> const& M1); /**< Operator overloading for subtraction */ 	
	fdcl::FFTSO3_matrix<ScalarType> operator-(fdcl::FFTSO3_matrix<double> const& M2); /**< Operator overloading for subtraction */ 		
    
    fdcl::FFTSO3_matrix<complex<double>> operator*(const complex<double>& c);  /**< Operator overloading for multiplying a complex number to each element */	
    fdcl::FFTSO3_matrix<ScalarType> operator*(const double& c);  /**< Operator overloading for multiplying a real number to each element */ 	

private:
	void assert_index(int l);
	void assert_index(int l, int m, int n);
};

namespace fdcl
{
    typedef fdcl::FFTSO3_matrix<double> FFTSO3_matrix_real; /**< Template implementation of fdcl::FFTSO3_matrix for real elements */ 
    typedef fdcl::FFTSO3_matrix<complex<double>> FFTSO3_matrix_complex; /**< Template implemmentation of fdcl::FFTSO3_matrix for complex elements */
}

#endif
