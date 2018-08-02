#ifndef _FDCL_FFTSO3_COMPLEX_HPP
#define _FDCL_FFTSO3_COMPLEX_HPP

#include <iostream>
#include <vector>
#include <math.h> // pow
#include <Eigen/Dense>

#include "fdcl_tictoc.hpp"
#include "fdcl_Clebsch_Gordon_matrix.hpp"
#include "misc_matrix_func.h"

using namespace std;
using namespace Eigen;

typedef Eigen::Matrix<double, 3, 3> Matrix3;

class fdcl_FFTSO3_complex
{
    public:
        std::vector<fdcl_FFTSO3_matrix_real> d_beta;
        fdcl_FFTSO3_matrix_complex D, F, F0, u;
        int B, l_max;
        std::vector<double> weight;
        fdcl_Clebsch_Gordon_matrix C;

        fdcl_FFTSO3_complex(){};
        fdcl_FFTSO3_complex(int l_max);
        ~fdcl_FFTSO3_complex(){};

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
        complex<double> inverse_transform(double alpha, double beta, double gamma);
        complex<double> inverse_transform(Matrix3);

        std::vector<fdcl_FFTSO3_matrix_complex> deriv_D();
        std::vector<double> character(double beta);

        // test
        void check_weight();
        void check_wigner_d();
        void check_deriv_D();
        void check_forward_transform();
        void check_Clebsch_Gordon();

    private:
        fdcl_FFTSO3_matrix_complex F4check;
        double delta(int ,int );
        double beta_k(int k);
        double alpha_j(int j);
        double gamma_j(int j);
        Eigen::VectorXd Legendre_poly(double x, int n);
        int signum(int );
        std::vector<double> compute_weight();

        static complex<double> f_4_check_forward_transform(double alpha, double beta, double gamma);

        fdcl_FFTSO3_matrix_complex forward_transform_0(std::function <complex<double>(double, double, double)>);

};

#endif
