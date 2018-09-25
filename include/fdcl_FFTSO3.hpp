#ifndef _FDCL_FFTSO3_HPP
#define _FDCL_FFTSO3_HPP

#include <iostream>
#include <vector>
#include <math.h> // pow
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

#include "fdcl_tictoc.hpp"
#include "fdcl_FFTSO3_matrix.hpp"
#include "fdcl_Clebsch_Gordon.hpp"
#include "fdcl_omp_thread.hpp"
#include "misc_matrix_func.h"

namespace fdcl
{
    class FFTSO3_complex;
    class FFTSO3_real;
}

/** \brief Complex Fast Fourier Transform on SO(3)
 *
 * This class provides various tools for complex harmonic analysis on the speical orthogonal group, includeing fast Forward transform, inverse transform, wigner-D function, wigner-d function, and Clebsch-Gordon coefficient. 
 *
 */
class fdcl::FFTSO3_complex
{
    public:
        /** \brief The maximum order of the Fourier tranform that is applied to all of the member functions. 
         *
         * For example, the forward trasnform is computed upto \f$l\leq l_{\max}\f$, i.e., \f$ F^{l}_{m,n} =0 \mbox{ for } l> l_{\max} \f$. Thhis can be set by the constructor, or fdcl::FFTSO3_complex::init(int l_max).
         */
        int l_max; 
        /** \brief fdcl::Clebsch_Gordon_complex class object
         *
         * Class instance for Clebsch-Gordon coefficients, which are computed by fdcl::Clebsch_Gordon_complex::compute(int l1, int l2). Then, the coefficient \f$ C^{l,m}_{l_1,m_1,l_2,m_2}\f$ is accessed by \c C(l,m,l1,m1,l2,m2).
         */
        fdcl::Clebsch_Gordon_complex C; 

        FFTSO3_complex(){};
        /** Constructor with \f$ l_{\max} \f$
         */
        FFTSO3_complex(int l_max);
        ~FFTSO3_complex(){};
        /** Initialize the class with the maximum order specificed by \c l_max.
         */
        void init(int l_max);

        /** wigner-d matrix
         *
         * Compute and returns a fdcl::FFTSO3_matrix_real object for real-valued wigner-d matrix, namely \f$ d^l_{m,n}(\beta) \f$. 
         * Each element can be accesed by the index  \c (l,m,n) of the returned object, or the \f$ (2l+1)\tims (2l+1)\f$ matrix for the \f$l\f$-th order is accessed by \c [l].
         */
        fdcl::FFTSO3_matrix_real wigner_d(double beta);

        /** wigner-D matrix for a given rotation matrix
         *
         * Compute and returns a fdcl::FFTSO3_matrix_complex object for complex-valued wigner-D matrix, namely \f$ D^l_{m,n}(R) \f$. 
         * Each element can be accesed by the index  \c (l,m,n) of the returned object, or the \f$ (2l+1)\times (2l+1)\f$ matrix for the \f$l\f$-th order is accessed by \c [l].
         */
        fdcl::FFTSO3_matrix_complex wigner_D(Eigen::Matrix3d);	
        /** wigner-D matrix for a given 3-2-3 Euler angles
         *
         * Compute and return a fdcl::FFTSO3_matrix_complex object for complex-valued wigner-D matrix, namely \f$ D^l_{m,n}(\alpha,\beta,\gamma) \f$. 
         * Each element can be accesed by the index  \c (l,m,n) of the returned object, or the \f$ (2l+1)\times (2l+1)\f$ matrix for the \f$l\f$-th order is accessed by \c [l].
         */
        fdcl::FFTSO3_matrix_complex wigner_D(double alpha, double beta, double gamma);

        /** Forward transform for \f$ f(\alpha,\beta,\gamma)\f$
         *
         * Compute and return a fdcl::FFTSO3_matrix_complex object for complex Fourier transform of \f$ f(\alpha,\beta,\gamma)):\mathrm{SO(3)}\rightarrow \mathbb{C} \f$, namely \f$ F^l_{m,n} \f$.
         * The function is defined in terms of 3-2-3 Euler angles. 
         * Each element can be accesed by the index  \c (l,m,n) of the returned object, or the \f$ (2l+1)\times (2l+1)\f$ matrix for the \f$l\f$-th order is accessed by \c [l].
         */
        fdcl::FFTSO3_matrix_complex forward_transform(std::function <complex<double>(double, double, double)>);

        /** Forward transform for \f$ f(R)\f$
         *
         * Compute and return a fdcl::FFTSO3_matrix_complex object for complex Fourier transform of \f$ f(R):\mathrm{SO(3)}\rightarrow \mathbb{C} \f$, namely \f$ F^l_{m,n} \f$.
         * The function is defined in terms of a rotation matrix. 
         * Each element can be accesed by the index  \c (l,m,n) of the returned object, or the \f$ (2l+1)\times (2l+1)\f$ matrix for the \f$l\f$-th order is accessed by \c [l].
         */
        fdcl::FFTSO3_matrix_complex forward_transform(std::function <complex<double>(Eigen::Matrix3d)>);

        /** Inverse trasnform for \f$(\alpha,\beta,\gamma)\f$
         *
         * Compute and return a complex number for the Fourer parameter evaluated at the 3-2-3 Euler angles \f$ (\alpha,\beta,\gamma) \f$.
         */
        complex<double> inverse_transform(fdcl::FFTSO3_matrix_complex, double alpha, double beta, double gamma);

        /** Inverse trasnform for \f$R\f$
         *
         * Compute and return a complex number for the Fourer parameter evaluated at the rotation matrix \f$ R \f$.
         */
        complex<double> inverse_transform(fdcl::FFTSO3_matrix_complex, Eigen::Matrix3d);

        /** Derivaties of wigner-D matrix
         *
         * Compute and return the derivatives of the wigner-D matrix at the identity along the canonial basis.
         * More specifically, let the output object be \c u. 
         * Then \c u[1], \c u[2], and \c u[3] correspond to \f$ u^l(e_1) \f$, \f$u^l(e_2)\f$, and \f$u^l(e_3)\f$, respectively. 
         * Note \c u[0] is set to zero.
         */
        std::vector<fdcl::FFTSO3_matrix_complex> deriv_wigner_D();

        /** Character 
         *
         * Compute and return the character \f$\chi(\beta)\f$. 
         */
        std::vector<double> character(double beta);

        /**  When set to \c true, the detailed results of check functions are printed. */
        bool check_verbose=false;         

        /** Execute all of check functions */
        void check_all();

        /** Check the properties of the weight used in fast Fourier transform 
         *
         * Verify the following two properties of \f$ w_k\f$. 
         * \f{align*}{
            \sum_k w_k d^l_{0,0}(\beta_k) 4B^2 & = \delta_{0,l},\\
            \sum_{j_1,k,j_2} w_k D(\alpha_{j_1}, \beta_k, \gamma_{j_2}) & = \delta_{l,0}\delta_{m,0}\delta_{n,0}
            \f}
        */
        double check_weight();

        /** Check the properties of the wigner-d function, computed by fdcl::FFTSO3_complex::wigner_d()
         *
         * Verify the orthogonalilty, namely \f$ d^l(\beta)(d^l(\beta))^T = I \f$, and the difference explicit expressions is at the order of machine precision.
         */
        double check_wigner_d();

        /** Chech the derivatives of the wigner-D function, computed by fdcl::FFTSO3_complex::deriv_wigner_D()
         *
         * Verify that the numerical derivative of the wigner-D function obtained by a finite-differece rule is consistent with fdcl::FFTSO3_complex::deriv_wigner_D()
         */
        double check_deriv_wigner_D();

        /** Check the forward transform and the inverse transform. 
         *
         * Verify that the composition of the inverse transform and the foward transform yield the identify map in the space of Fourier parameters. 
         * More specifically, a function is defined as the form of the inverse trasnform with random Fourier parameters, and check if its Fourier transform yields the same Fourier parameters. 
         */
        double check_transform();
        
        /** Check the Clebsh-Gordon coefficients, computed by fdcl::Clebsch_Gordon_complex::compute()
         *
         * Verify the followign property:
         * \f{equation*}{
         D^{l_1}_{m_1,n_1} (R) D^{l_2}_{m_2,n_2}(R) = \sum_{l=\underline{l}}^{l_1+l_2} C^{l,m_1+m_2}_{l_1,m_1,l_2,m_2} C^{l,n_1+n_2}_{l_1,n_1,l_2,n_2} D^{l}_{m_1+m_2,n_1+n_2}(R).
         \f}
         */
        double check_Clebsch_Gordon();


	protected:
        std::vector<double> weight;
        std::vector<fdcl::FFTSO3_matrix_real> d_beta;
        int B;  
        fdcl::FFTSO3_matrix_complex D, F, u;
        double beta_k(int k);
        double alpha_j(int j);
        double gamma_j(int j);
        std::vector<double> compute_weight();
        fdcl::FFTSO3_matrix_real wigner_d(double beta, int L);	
        fdcl::FFTSO3_matrix_real wigner_d_explicit(double beta);
        fdcl::FFTSO3_matrix_complex wigner_D(double alpha, double beta, double gamma, int L);
        fdcl::FFTSO3_matrix_complex forward_transform(std::function <complex<double>(double, double, double)>, bool is_real);

    private:
        fdcl::FFTSO3_matrix_complex F_4_check;
        complex<double> f_4_check_transform(double alpha, double beta, double gamma);
		fdcl::FFTSO3_matrix_complex forward_transform_0(std::function <complex<double>(double, double, double)>);
        fdcl::FFTSO3_matrix_complex forward_transform_1(std::function <complex<double>(double, double, double)>);
};

/** \brief Real Fast Fourier Transform on SO(3)
 *
 * This class provides various tools for real harmonic analysis on the speical orthogonal group, includeing fast Forward transform, inverse transform, wigner-D function, wigner-d function, and Clebsch-Gordon coefficient. 
 *
 */
class fdcl::FFTSO3_real : public fdcl::FFTSO3_complex
{
    public:
        /** \brief fdcl::Clebsch_Gordon_real class object
         *
         * Class instance for Clebsch-Gordon coefficients, which are computed by fdcl::Clebsch_Gordon_real::compute(int l1, int l2). Then, the coefficient \f$ c^{l,m}_{l_1,m_1,l_2,m_2}\f$ is accessed by \c c(l,m,l1,m1,l2,m2).
         */
        fdcl::Clebsch_Gordon_real c;

        FFTSO3_real() {};
        /** Constructor with \f$ l_{\max} \f$
         */
        FFTSO3_real(int l_max);
        ~FFTSO3_real(){};

        /** Initialize the class with the maximum order specificed by \c l_max.
        */ 
        void init(int l_max);

        /** real harmonics on SO(3) for given 3-2-3 Euler angles
         *
         * Compute and returns a fdcl::FFTSO3_matrix_real object for real harmonics on SO(3), namely \f$ U^l_{m,n}(R) \f$. 
         * Each element can be accesed by the index  \c (l,m,n) of the returned object, or the \f$ (2l+1)\times (2l+1)\f$ matrix for the \f$l\f$-th order is accessed by \c [l].
         */
        fdcl::FFTSO3_matrix_real real_harmonics(double alpha, double beta, double gamma);

        /** real harmonics on SO(3) for a given rotation matrix
         *
         * Compute and returns a fdcl::FFTSO3_matrix_real object for real harmonics on SO(3), namely \f$ U^l_{m,n}(\alpha,\beta,\gamma) \f$. 
         * Each element can be accesed by the index  \c (l,m,n) of the returned object, or the \f$ (2l+1)\times (2l+1)\f$ matrix for the \f$l\f$-th order is accessed by \c [l].
         */
        fdcl::FFTSO3_matrix_real real_harmonics(Eigen::Matrix3d);

        /** Derivaties of wigner-D matrix
         *
         * Compute and return the derivatives of the real harmonics at the identity along the canonial basis.
         * More specifically, let the output object be \c u. 
         * Then \c u[1], \c u[2], and \c u[3] correspond to \f$ u^l(e_1) \f$, \f$u^l(e_2)\f$, and \f$u^l(e_3)\f$, respectively. 
         * Note \c u[0] is set to zero.
         */
        std::vector<fdcl::FFTSO3_matrix_real> deriv_real_harmonics();

        /** Forward transform for \f$ f(\alpha,\beta,\gamma)\f$
         *
         * Compute and return a fdcl::FFTSO3_matrix_real object for real Fourier transform of \f$ f(\alpha,\beta,\gamma)):\mathrm{SO(3)}\rightarrow \mathbb{C} \f$, namely \f$ F^l_{m,n} \f$.
         * The function is defined in terms of 3-2-3 Euler angles. 
         * Each element can be accesed by the index  \c (l,m,n) of the returned object, or the \f$ (2l+1)\times (2l+1)\f$ matrix for the \f$l\f$-th order is accessed by \c [l].
         */
        fdcl::FFTSO3_matrix_real forward_transform(std::function <double(double, double, double)>);

        /** Forward transform for \f$ f(R)\f$
         *
         * Compute and return a fdcl::FFTSO3_matrix_real object for real Fourier transform of \f$ f(R):\mathrm{SO(3)}\rightarrow \mathbb{C} \f$, namely \f$ F^l_{m,n} \f$.
         * The function is defined in terms of a rotation matrix. 
         * Each element can be accesed by the index  \c (l,m,n) of the returned object, or the \f$ (2l+1)\times (2l+1)\f$ matrix for the \f$l\f$-th order is accessed by \c [l].
         */
        fdcl::FFTSO3_matrix_real forward_transform(std::function <double(Eigen::Matrix3d)>);

        /** Inverse trasnform for \f$(\alpha,\beta,\gamma)\f$
         *
         * Compute and return a real number for the Fourer parameter evaluated at the 3-2-3 Euler angles \f$ (\alpha,\beta,\gamma) \f$.
         */
        double inverse_transform(fdcl::FFTSO3_matrix_real, double alpha, double beta, double gamma);

        /** Inverse trasnform for \f$R\f$
         *
         * Compute and return a real number for the Fourer parameter evaluated at the rotation matrix \f$ R \f$.
         */
        double inverse_transform(fdcl::FFTSO3_matrix_real, Eigen::Matrix3d);

        /** Execute all of check functions */
        void check_all();

        /** Check real harmonics \f$U^l_{m,n}(R) \f$
         *
         * Verify that the real harmonics computed by various methods are consistent with each other
         */
        double check_real_harmonics();

        /** Check the Clebsh-Gordon coefficients, computed by fdcl::Clebsch_Gordon_complex::compute()
         *
         * Verify the followign property:
         * \f{equation*}{
         U^{l_1}_{m_1,n_1} (R) U^{l_2}_{m_2,n_2}(R) =  \sum_{m\in M} \sum_{n\in N} \sum_{l=\max\{|l_1-l_2|,|m|,|n|\}}^{l_1+l_2} c^{l,m}_{l_1,m_1,l_2,m_2} \overline{c}^{l,n}_{l_1,n_1,l_2,n_2} U^{l}_{m,n}(R).
         \f}
         */
        double check_Clebsch_Gordon();

        /** Chech the derivatives of the real harmonics, computed by fdcl::FFTSO3_real::real_harmonics()
         *
         * Verify that the numerical derivative of the real harmonics obtained by a finite-differece rule is consistent with fdcl::FFTSO3_real::real_harmonics()
         */
        double check_deriv_real_harmonics();

        /** Check the forward transform and the inverse transform. 
         *
         * Verify that the composition of the inverse transform and the foward transform yield the identify map in the space of Fourier parameters. 
         * More specifically, a function is defined as the form of the inverse trasnform with random Fourier parameters, and check if its Fourier transform yields the same Fourier parameters. 
         */
        double check_transform();

    private:
        fdcl::FFTSO3_matrix_real U;
        fdcl::FFTSO3_matrix_real real_harmonics(double alpha, double beta, double gamma, int L);
        int index_fft(int , int);
		int signum(int );
        fdcl::FFTSO3_matrix_real F_4_check;
        double f_4_check_transform(double alpha, double beta, double gamma);

        fdcl::FFTSO3_matrix_real real_harmonics_3(double alpha, double beta, double gamma, int L); // alternative method with matrix formulation
        fdcl::FFTSO3_matrix_complex real_harmonics_2(double alpha, double beta, double gamma, int L); // alternative method with U = \bar C D C^T
        fdcl::FFTSO3_matrix_real real_harmonics_1(double alpha, double beta, double gamma, int L);// alternative formulation based on Phi_1 and Phi_2
        fdcl::FFTSO3_matrix_real real_harmonics_0(double alpha, double beta, double gamma, int L);// alternative formulation based on Psi

        fdcl::FFTSO3_matrix_real forward_transform_0(std::function <double(double, double, double)>);
        fdcl::FFTSO3_matrix_real forward_transform_1(std::function <double(double, double, double)>);
        fdcl::FFTSO3_matrix_real forward_transform_2(std::function <double(double, double, double)>);

        std::vector<double> compute_Phi(int m, int n, double alpha, double gamma);	
        fdcl::FFTSO3_matrix_real compute_Psi(double beta, int L);
        void compute_X(double gamma, int L, Eigen::MatrixXd& X);

        fdcl::FFTSO3_matrix_complex T;
        fdcl::FFTSO3_matrix_complex matrix2rsph(int L);
};

#endif
