#ifndef _FDCL_CLEBSCH_GORDON_HPP
#define _FDCL_CLEBSCH_GORDON_HPP

#include <iostream>
#include <vector>
#include <math.h> // pow
#include <array>
#include <Eigen/Dense>

#include "fdcl_FFTSO3_matrix.hpp"
#include "fdcl_tictoc.hpp"

namespace fdcl
{
    class Clebsch_Gordon_complex;
    class Clebsch_Gordon_real;
}

/** \brief Clebsh-Gordon Coefficient for Complex Harmonics on SO(3)
 *
 * This class formulates a matrix class for the Clebsh-Gordon coefficients for complex harmonics on SO(3)
 * The Clebsh-Gordon coefficients are indexed by six integers, i.e., \f$ C^{l,m}_{l_1,m_1,l_2,m_2}\f$.
 * This class overload the operator \c() such that the above element can be accessed by \c C(l,m,l1,m1,l2,n2).
 *
 * NOTE:
 *  - Prior to accessing the element, the member function \c compute(l1,l2) must be called with the specific value of \f$ (l_1,l_2) \f$.
 *  - The elements for other values of \f$ (l,m,m_1,m_2) \f$ with the same \f$ (l_1,l_2)\f$ as above can be accessed without need for calling \c compute(l1,l2)
 *  - However, to access the elements with another value of \f$ (l'_1,l'_2) \f$, the member function must be called again with \c compute(l1_prime, l2_prime)
 *
 */
class fdcl::Clebsch_Gordon_complex
{
    public:
        Clebsch_Gordon_complex(){};
        Clebsch_Gordon_complex(int l1, int l2); /**< Constructor with \c l1 and \c l2 */
        void init(int l1, int l2); /**< Initialize with \c l1 and \c l2 */
        ~Clebsch_Gordon_complex(){};
        int l1, l2;
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> C; /**< Buffer */

        double& operator()(int l, int m, int l1, int m1, int l2, int m2); /**< Operator overloading to access \f$ C^{l,m}_{l_1,m_1,l_2,m_2}\f$.*/ 

        void compute(int l1, int l2);  /**< Compute all of the coefficients with the given \f$l_1,l_2\f$. Coefficients for other values of \f$l_1,l_2\f$ are NOT computed. */
    protected:
        int row(int l, int m, int l1, int m1, int l2, int m2);
        int col(int l, int m, int l1, int m1, int l2, int m2);
        void assert_index(int l, int m, int l1, int m1, int l2, int m2);
        void compute_sub(int l, int m, int l1, int l2);
        fdcl::FFTSO3_matrix_complex matrix2rsph(int );
};

/** \brief Clebsh-Gordon Coefficient for Real Harmonics on SO(3)
 *
 * This class formulates a matrix class for the Clebsh-Gordon coefficients for real harmonics on SO(3)
 * The Clebsh-Gordon coefficients are indexed by six integers, i.e., \f$ c^{l,m}_{l_1,m_1,l_2,m_2}\f$.
 * This class overload the operator \c() such that the above element can be accessed by \c c(l,m,l1,m1,l2,n2)
 *
 * NOTE:
 *  - Prior to accessing the element, the member function \c compute(l1,l2) must be called with the specific value of \f$ (l_1,l_2) \f$.
 *  - The elements for other values of \f$ (l,m,m_1,m_2) \f$ with the same \f$ (l_1,l_2)\f$ as above can be accessed without need for calling \c compute(l1,l2)
 *  - However, to access the elements with another value of \f$ (l'_1,l'_2) \f$, the member function must be called again with \c compute(l1_prime, l2_prime)
 *
 */
class fdcl::Clebsch_Gordon_real : public fdcl::Clebsch_Gordon_complex
{
    public:
        Clebsch_Gordon_real(){};
        Clebsch_Gordon_real(int l1, int l2); /**< Constructor with \c l1 and \c l2 */
        void init(int l1, int l2); /**< Initialize with \c l1 and \c l2 */
        ~Clebsch_Gordon_real(){};

        Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> c; /**< Buffer */
        complex<double>& operator()(int l, int m, int l1, int m1, int l2, int m2);  /**< Operator overloading to access \f$ C^{l,m}_{l_1,m_1,l_2,m_2}\f$.*/ 
        void compute(int l1, int l2);  /**< Compute all of the coefficients with the given \f$l_1,l_2\f$. Coefficients for other values of \f$l_1,l_2\f$ are NOT computed. */
        void print();
    private:
        Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> X;
        void compute_0(int l1, int l2);
        void compute_sub_01(int l1, int m1, int l2, int m2, complex<double> ratio0, complex<double> ratio1);
        void compute_sub_23(int l1, int m1, int l2, int m2, complex<double> ratio2, complex<double> ratio3);

};

#endif
