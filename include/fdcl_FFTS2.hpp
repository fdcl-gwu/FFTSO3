#ifndef _FDCL_FFTS2_HPP
#define _FDCL_FFTS2_HPP

#include <iostream>
#include <vector>
#include <math.h> // pow
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <omp.h>

#include "fdcl_tictoc.hpp"
#include "fdcl_FFTSO3_matrix.hpp"
#include "fdcl_FFTS2_matrix.hpp"
#include "misc_matrix_func.h"
#include "fdcl_omp_thread.hpp"

namespace fdcl
{
    class FFTS2_complex;
    class FFTS2_real;
}

/** \brief Complex Fast Fourier Transform on Sphere
 *
 * This class provides various tools for complex harmonic analysis on the two sphere, includeing fast Forward transform, inverse transform.
 *
 */
class fdcl::FFTS2_complex
{
    public:
        /** The maximum order of the Fourier tranform that is applied to the foward/inverse transform.
         */
        int l_max;
        FFTS2_complex(){};

        /** Constructor with \f$ l_{\max} \f$
         */
        FFTS2_complex(int l_max);
        ~FFTS2_complex(){};
        /** Initialize the class with the maximum order specificed by \c l_max.
         */
        void init(int l_max);

        /** Spherical harmonics 
         *
         * Compute and return a fdcl::FFTS2_matrix_complex object for complex-valued spherical harmonics, namely \f$ Y^l_{m}(\theta,\phi) \f$ upto \f$ l\leq L\f$. 
         * \f[
         * Y^l_m(\theta,\phi) = e^{im\phi} \sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}}  P^m_l(\cos\theta) = e^{im\phi}\Theta^l_m(\theta).
         \f]
         * Each element can be accesed by the index  \c (l,m) of the returned object.
         */
        fdcl::FFTS2_matrix_complex spherical_harmonics(double theta, double phi, int L);

        /** Forward transform for \f$ f(\theta,\phi)\f$
         *
         * Compute and return a fdcl::FFTS2_matrix_complex object for complex Fourier transform of \f$ f(\theta,\phi):\mathrm{S^2}\rightarrow \mathbb{C} \f$, namely \f$ F^l_{m} \f$.
         * The function is defined in terms of the polar coordinates \f$(\theta,\phi)\in[0,\pi]\times[0,2\pi)\f$, i.e., \f$x\in\mathrm{S}^2\f$ is parameterized by \f$x=(\cos\phi\sin\theta, \sin\phi\sin\theta, \cos\theta)\f$. 
         */
        fdcl::FFTS2_matrix_complex forward_transform(std::function <complex<double>(double, double)>);

        /** Inverse trasnform for \f$R\f$
         *
         * Compute and return a complex number for the Fourer parameter evaluated at the coordinate \f$(\theta,\phi)\f$.
         */
        complex<double> inverse_transform(fdcl::FFTS2_matrix_complex F, double theta, double phi);

        /**  When set to \c true, the detailed results of check functions are printed. */
        bool check_verbose=false;

        /** Execute all of check functions */
        void check_all();
        /** Check the properties of the weight used in fast Fourier transform 
       */
        double check_weight();

        /** Check the forward transform and the inverse transform. 
         *
         * Verify that the composition of the inverse transform and the foward transform yield the identify map in the space of Fourier parameters. 
         * More specifically, a function is defined as the form of the inverse trasnform with random Fourier parameters, and check if its Fourier transform yields the same Fourier parameters. 
         */
        double check_transform();

     protected:
        int B;
        fdcl::FFTS2_matrix_complex Y, F;
        std::vector<double> weight;

        complex<double> inverse_transform(double theta, double phi);
        fdcl::FFTS2_matrix_complex forward_transform(std::function <complex<double>(double, double)>, bool );
        fdcl::FFTS2_matrix_real nor_assoc_Legendre_poly(double cos_beta, int L);
        std::vector<double> compute_weight();
        double theta_k(int k);
        double phi_j(int j);

        fdcl::FFTS2_matrix_real nP;

    private:
        fdcl::FFTS2_matrix_complex F_4_check;
        complex<double> f_4_check_transform(double theta, double phi);
};

/** \brief Real Fast Fourier Transform on Sphere
 *
 * This class provides various tools for real harmonic analysis on the two sphere, includeing fast Forward transform, inverse transform.
 *
 */
class fdcl::FFTS2_real : public fdcl::FFTS2_complex
{
    public:
        FFTS2_real(){};

        /** Constructor with \f$ l_{\max} \f$
         */
        FFTS2_real(int l_max);

        ~FFTS2_real(){};
        /** Initialize the class with the maximum order specificed by \c l_max.
         */
        void init(int l_max);

        /** Real spherical harmonics 
         *
         * Compute and return a fdcl::FFTS2_matrix_real object for real-valued spherical harmonics, namely \f$ S^l_{m}(\theta,\phi) \f$ upto \f$ l\leq L\f$. 
         * This function follows the convention of real harmonics presented in 
         *    - M. Blanco, M. Fl´orez, and M. Bermejo, “Evaluation of the rotation matrices in the basis of real spherical harmonics,” Journal of
Molecular Structure, vol. 419, pp. 19–27, 1997.
         * Each element can be accesed by the index  \c (l,m) of the returned object.
         */
        fdcl::FFTS2_matrix_real spherical_harmonics(double theta, double phi, int L);

        /** Forward transform for \f$ f(\theta,\phi)\f$
         *
         * Compute and return a fdcl::FFTS2_matrix_real object for real Fourier transform of \f$ f(\theta,\phi):\mathrm{S^2}\rightarrow \mathbb{C} \f$, namely \f$ F^l_{m} \f$.
         * The function is defined in terms of the polar coordinates \f$(\theta,\phi)\in[0,\pi]\times[0,2\pi)\f$, i.e., \f$x\in\mathrm{S}^2\f$ is parameterized by \f$x=(\cos\phi\sin\theta, \sin\phi\sin\theta, \cos\theta)\f$. 
         */
        fdcl::FFTS2_matrix_real forward_transform(std::function <double(double, double)>);

        /** Inverse trasnform for \f$R\f$
         *
         * Compute and return a complex number for the Fourer parameter evaluated at the coordinate \f$(\theta,\phi)\f$.
         */
        double inverse_transform(fdcl::FFTS2_matrix_real F, double theta, double phi);
       
        /** Execute all of check functions */
        void check_all();

        /** Check the forward transform and the inverse transform. 
         *
         * Verify that the composition of the inverse transform and the foward transform yield the identify map in the space of Fourier parameters. 
         * More specifically, a function is defined as the form of the inverse trasnform with random Fourier parameters, and check if its Fourier transform yields the same Fourier parameters. 
         */
        double check_transform();

    private:
        fdcl::FFTSO3_matrix_complex T;
        fdcl::FFTSO3_matrix_complex matrix2rsph(int L);
        fdcl::FFTS2_matrix_real y, F;
        fdcl::FFTS2_matrix_real F_4_check;
        double f_4_check_transform(double theta, double phi);
};

#endif
