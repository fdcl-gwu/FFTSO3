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
#include "misc_matrix_func.h"

using namespace std;
using namespace Eigen;

class fdcl_FFTSO3_complex
{
    public:
        std::vector<fdcl_FFTSO3_matrix_real> d_beta;
        fdcl_FFTSO3_matrix_complex D, F, u;
        int B, l_max;
        std::vector<double> weight;
        fdcl_Clebsch_Gordon_complex C;

        fdcl_FFTSO3_complex(){};
        fdcl_FFTSO3_complex(int l_max);
        ~fdcl_FFTSO3_complex(){};
        void init(int l_max);

        fdcl_FFTSO3_matrix_real wigner_d(double beta, int L);	
        fdcl_FFTSO3_matrix_real wigner_d(double beta);
        fdcl_FFTSO3_matrix_real wigner_d_explicit(double beta);

        // complex transform
        fdcl_FFTSO3_matrix_complex wigner_D(Matrix3);	
        fdcl_FFTSO3_matrix_complex wigner_D(double alpha, double beta, double gamma);
        fdcl_FFTSO3_matrix_complex wigner_D(double alpha, double beta, double gamma, int L);

        fdcl_FFTSO3_matrix_complex forward_transform(std::function <complex<double>(double, double, double)>);
        fdcl_FFTSO3_matrix_complex forward_transform(std::function <complex<double>(Matrix3)>);

        complex<double> inverse_transform(fdcl_FFTSO3_matrix_complex, double alpha, double beta, double gamma);
        complex<double> inverse_transform(fdcl_FFTSO3_matrix_complex, Matrix3);

        std::vector<fdcl_FFTSO3_matrix_complex> deriv_wigner_D();
        std::vector<double> character(double beta);

        // test
        void check_weight();
        void check_wigner_d();
        void check_deriv_wigner_D();
        void check_transform();
        void check_Clebsch_Gordon();

	protected:
        double beta_k(int k);
        double alpha_j(int j);
        double gamma_j(int j);
        std::vector<double> compute_weight();

    private:
        fdcl_FFTSO3_matrix_complex F_4_check;
        complex<double> f_4_check_transform(double alpha, double beta, double gamma);
		fdcl_FFTSO3_matrix_complex forward_transform_0(std::function <complex<double>(double, double, double)>);
        fdcl_FFTSO3_matrix_complex forward_transform_1(std::function <complex<double>(double, double, double)>);
};

class fdcl_FFTSO3_real : public fdcl_FFTSO3_complex
{
    public:
        fdcl_Clebsch_Gordon_real c;
        fdcl_FFTSO3_matrix_real U;

        fdcl_FFTSO3_real() {};
        fdcl_FFTSO3_real(int l_max);
        ~fdcl_FFTSO3_real(){};
        void init(int l_max);

        // real transform
        fdcl_FFTSO3_matrix_real wigner_D_real(double alpha, double beta, double gamma, int L);
        fdcl_FFTSO3_matrix_real wigner_D_real(double alpha, double beta, double gamma);
        fdcl_FFTSO3_matrix_real wigner_D_real(Matrix3);

        std::vector<fdcl_FFTSO3_matrix_real> deriv_wigner_D_real();
        fdcl_FFTSO3_matrix_real forward_transform(std::function <double(double, double, double)>);
        fdcl_FFTSO3_matrix_real forward_transform(std::function <double(Matrix3)>);

        double inverse_transform(fdcl_FFTSO3_matrix_real, double alpha, double beta, double gamma);
        double fast_inverse_transform(fdcl_FFTSO3_matrix_real, double alpha, double beta, double gamma);
        double inverse_transform(fdcl_FFTSO3_matrix_real, Matrix3);

        void check_wigner_D_real();
        void check_Clebsch_Gordon();
        void check_deriv_wigner_D_real();
        void check_transform();

    private:
		int signum(int );
        fdcl_FFTSO3_matrix_real F_4_check;
        double f_4_check_transform(double alpha, double beta, double gamma);

        fdcl_FFTSO3_matrix_complex wigner_D_real_2(double alpha, double beta, double gamma, int L); // alternative method with U = \bar C D C^T
        fdcl_FFTSO3_matrix_real wigner_D_real_1(double alpha, double beta, double gamma, int L);// alternative formulation based on Phi_1 and Phi_2
        fdcl_FFTSO3_matrix_real wigner_D_real_0(double alpha, double beta, double gamma, int L);// alternative formulation based on Theta/Psi

        fdcl_FFTSO3_matrix_real forward_transform_0(std::function <double(double, double, double)>);
        fdcl_FFTSO3_matrix_real forward_transform_1(std::function <double(double, double, double)>);

        std::vector<double> compute_Phi(int m, int n, double alpha, double gamma);	
        std::vector<fdcl_FFTSO3_matrix_real> compute_Theta_Psi(double beta, int L);

        fdcl_FFTSO3_matrix_complex T;
        fdcl_FFTSO3_matrix_complex matrix2rsph(int L);
};

#endif
