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
	fdcl_FFTSO3_matrix_complex wigner_D_real(double alpha, double beta, double gamma, int L);
	fdcl_FFTSO3_matrix_complex wigner_D_real(double alpha, double beta, double gamma);
	fdcl_FFTSO3_matrix_complex wigner_D_real(Matrix3);
    fdcl_FFTSO3_matrix_real wigner_D_real_direct(double alpha, double beta, double gamma, int L);
	fdcl_FFTSO3_matrix_real wigner_D_real_direct(double alpha, double beta, double gamma);
	fdcl_FFTSO3_matrix_real wigner_D_real_direct(Matrix3);


    fdcl_FFTSO3_matrix_complex matrix2rsph(int L);

	complex<double> inverse_transform(fdcl_FFTSO3_matrix_complex, double alpha, double beta, double gamma);
	complex<double> inverse_transform(fdcl_FFTSO3_matrix_complex, Matrix3);
	complex<double> inverse_transform(double alpha, double beta, double gamma);
	complex<double> inverse_transform(Matrix3);
	
	fdcl_FFTSO3_matrix_complex forward_transform_0();
	fdcl_FFTSO3_matrix_complex forward_transform_1();

	std::vector<double> compute_weight();
	std::vector<fdcl_FFTSO3_matrix_complex> deriv_D();
	std::vector<double> character(double beta);

	void check_weight();
	void check_wigner_d();
	void check_deriv_D();

	complex<double> f(double alpha, double beta, double gamma);

private:
	double delta(int ,int );
	double beta_k(int k);
	double alpha_j(int j);
	double gamma_j(int j);
	Eigen::VectorXd Legendre_poly(double x, int n);
	
	
};

#endif
