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
class fdcl::FFTSO3_complex
{
    public:
        std::vector<fdcl::FFTSO3_matrix_real> d_beta;
        fdcl::FFTSO3_matrix_complex D, F, u;
        int B, l_max;
        std::vector<double> weight;
        fdcl::Clebsch_Gordon_complex C;

        FFTSO3_complex(){};
        FFTSO3_complex(int l_max);
        ~FFTSO3_complex(){};
        void init(int l_max);

        fdcl::FFTSO3_matrix_real wigner_d(double beta, int L);	
        fdcl::FFTSO3_matrix_real wigner_d(double beta);
        fdcl::FFTSO3_matrix_real wigner_d_explicit(double beta);

        // complex transform
        fdcl::FFTSO3_matrix_complex wigner_D(Eigen::Matrix3d);	
        fdcl::FFTSO3_matrix_complex wigner_D(double alpha, double beta, double gamma);
        fdcl::FFTSO3_matrix_complex wigner_D(double alpha, double beta, double gamma, int L);

        fdcl::FFTSO3_matrix_complex forward_transform(std::function <complex<double>(double, double, double)>, bool is_real);
        fdcl::FFTSO3_matrix_complex forward_transform(std::function <complex<double>(double, double, double)>);
        fdcl::FFTSO3_matrix_complex forward_transform(std::function <complex<double>(Eigen::Matrix3d)>);

        complex<double> inverse_transform(fdcl::FFTSO3_matrix_complex, double alpha, double beta, double gamma);
        complex<double> inverse_transform(fdcl::FFTSO3_matrix_complex, Eigen::Matrix3d);

        std::vector<fdcl::FFTSO3_matrix_complex> deriv_wigner_D();
        std::vector<double> character(double beta);

        // test
        bool check_verbose=false;
        void check_all();
        double check_weight();
        double check_wigner_d();
        double check_deriv_wigner_D();
        double check_transform();
        double check_Clebsch_Gordon();

	protected:
        double beta_k(int k);
        double alpha_j(int j);
        double gamma_j(int j);
        std::vector<double> compute_weight();

    private:
        fdcl::FFTSO3_matrix_complex F_4_check;
        complex<double> f_4_check_transform(double alpha, double beta, double gamma);
		fdcl::FFTSO3_matrix_complex forward_transform_0(std::function <complex<double>(double, double, double)>);
        fdcl::FFTSO3_matrix_complex forward_transform_1(std::function <complex<double>(double, double, double)>);
};

class fdcl::FFTSO3_real : public fdcl::FFTSO3_complex
{
    public:
        fdcl::Clebsch_Gordon_real c;
        fdcl::FFTSO3_matrix_real U;

        FFTSO3_real() {};
        FFTSO3_real(int l_max);
        ~FFTSO3_real(){};
        void init(int l_max);

        // real transform
        fdcl::FFTSO3_matrix_real real_harmonics(double alpha, double beta, double gamma, int L);
        fdcl::FFTSO3_matrix_real real_harmonics(double alpha, double beta, double gamma);
        fdcl::FFTSO3_matrix_real real_harmonics(Eigen::Matrix3d);

        std::vector<fdcl::FFTSO3_matrix_real> deriv_real_harmonics();
        fdcl::FFTSO3_matrix_real forward_transform(std::function <double(double, double, double)>);
        fdcl::FFTSO3_matrix_real forward_transform(std::function <double(Eigen::Matrix3d)>);

        double inverse_transform(fdcl::FFTSO3_matrix_real, double alpha, double beta, double gamma);
        double fast_inverse_transform(fdcl::FFTSO3_matrix_real, double alpha, double beta, double gamma);
        double inverse_transform(fdcl::FFTSO3_matrix_real, Eigen::Matrix3d);

        void check_all();
        double check_real_harmonics();
        double check_Clebsch_Gordon();
        double check_deriv_real_harmonics();
        double check_transform();

    private:
		int signum(int );
        fdcl::FFTSO3_matrix_real F_4_check;
        double f_4_check_transform(double alpha, double beta, double gamma);

        fdcl::FFTSO3_matrix_complex real_harmonics_2(double alpha, double beta, double gamma, int L); // alternative method with U = \bar C D C^T
        fdcl::FFTSO3_matrix_real real_harmonics_1(double alpha, double beta, double gamma, int L);// alternative formulation based on Phi_1 and Phi_2
        fdcl::FFTSO3_matrix_real real_harmonics_0(double alpha, double beta, double gamma, int L);// alternative formulation based on Theta/Psi

        fdcl::FFTSO3_matrix_real forward_transform_0(std::function <double(double, double, double)>);
        fdcl::FFTSO3_matrix_real forward_transform_1(std::function <double(double, double, double)>);

        std::vector<double> compute_Phi(int m, int n, double alpha, double gamma);	
        fdcl::FFTSO3_matrix_real compute_Psi(double beta, int L);

        fdcl::FFTSO3_matrix_complex T;
        fdcl::FFTSO3_matrix_complex matrix2rsph(int L);
};

#endif
