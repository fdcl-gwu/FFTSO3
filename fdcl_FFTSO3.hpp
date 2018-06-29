#ifndef _FDCL_FFTSO3_HPP
#define _FDCL_FFTSO3_HPP

#include <iostream>
#include <vector>
#include <math.h> // pow
#include <Eigen/Dense>
#include "misc_matrix_func.h"

using namespace std;
using namespace Eigen;

typedef Eigen::Matrix<double, 3, 3> Matrix3;

class fdcl_FFTSO3
{
public:
	std::vector<fdcl_FFTSO3_matrix_real> d_beta;
	fdcl_FFTSO3_matrix_complex D, u, F, F0;
	int B, l_max;
	std::vector<double> weight;
	
	fdcl_FFTSO3(){};
	fdcl_FFTSO3(int l_max);
	~fdcl_FFTSO3(){};

	fdcl_FFTSO3_matrix_real wigner_d(double beta, int L);	
	fdcl_FFTSO3_matrix_real wigner_d(double beta);
	fdcl_FFTSO3_matrix_real wigner_d_explicit(double beta);

	fdcl_FFTSO3_matrix_complex wigner_D(Matrix3);	
	fdcl_FFTSO3_matrix_complex wigner_D(double alpha, double beta, double gamma);
	fdcl_FFTSO3_matrix_complex wigner_D(double alpha, double beta, double gamma, int L);

	fdcl_FFTSO3_matrix_real wigner_D_real(double alpha, double beta, double gamma, int L);
	fdcl_FFTSO3_matrix_real wigner_D_real(double alpha, double beta, double gamma);
	fdcl_FFTSO3_matrix_real wigner_D_real(Matrix3);
    fdcl_FFTSO3_matrix_complex wigner_D_real_converted(double alpha, double beta, double gamma, int L); // alternative method with U = \bar C D C^T
	fdcl_FFTSO3_matrix_real wigner_D_real_Phi(double alpha, double beta, double gamma, int L);// alternative formulation based on Phi_1 and Phi_2
    fdcl_FFTSO3_matrix_real wigner_D_real_0(double alpha, double beta, double gamma, int L);

	complex<double> inverse_transform(fdcl_FFTSO3_matrix_complex, double alpha, double beta, double gamma);
	complex<double> inverse_transform(fdcl_FFTSO3_matrix_complex, Matrix3);
	complex<double> inverse_transform(double alpha, double beta, double gamma);
	complex<double> inverse_transform(Matrix3);
	
	fdcl_FFTSO3_matrix_complex forward_transform_0();
	fdcl_FFTSO3_matrix_complex forward_transform_1();

    fdcl_FFTSO3_matrix_real forward_transform_real();
    fdcl_FFTSO3_matrix_real forward_transform_real_0();
    double inverse_transform(fdcl_FFTSO3_matrix_real, double alpha, double beta, double gamma);

	std::vector<fdcl_FFTSO3_matrix_complex> deriv_D();
	std::vector<double> character(double beta);

	void check_weight();
	void check_wigner_d();
	void check_deriv_D();

	complex<double> f(double alpha, double beta, double gamma);
	double f_real(double alpha, double beta, double gamma);

private:
	double delta(int ,int );
	double beta_k(int k);
	double alpha_j(int j);
	double gamma_j(int j);
	Eigen::VectorXd Legendre_poly(double x, int n);
    int signum(int );
    fdcl_FFTSO3_matrix_complex matrix2rsph(int L);
	std::vector<double> compute_weight();
    std::vector<double> compute_Phi(int m, int n, double alpha, double gamma);	
    std::vector<fdcl_FFTSO3_matrix_real> compute_Theta_Psi(double beta, int L);
};

#endif
