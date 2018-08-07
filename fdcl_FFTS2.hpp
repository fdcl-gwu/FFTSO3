#ifndef _FDCL_FFTS2_HPP
#define _FDCL_FFTS2_HPP

#include <iostream>
#include <vector>
#include <math.h> // pow
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

#include "fdcl_tictoc.hpp"
#include "fdcl_FFTSO3_matrix.hpp"
#include "fdcl_FFTS2_matrix.hpp"
#include "misc_matrix_func.h"

using namespace std;
using namespace Eigen;

class fdcl_FFTS2_complex
{
    public:
        int l_max, B;
        fdcl_FFTS2_matrix_complex Y;
        std::vector<double> weight;

        fdcl_FFTS2_complex(){};
        fdcl_FFTS2_complex(int l_max);
        ~fdcl_FFTS2_complex(){};
        void init(int l_max);

        fdcl_FFTS2_matrix_complex spherical_harmonics(double theta, double phi, int L);

        fdcl_FFTS2_matrix_complex forward_transform(std::function <complex<double>(double, double)>);
        fdcl_FFTS2_matrix_complex forward_transform(std::function <complex<double>(double, double)>, bool );
        complex<double> inverse_transform(fdcl_FFTS2_matrix_complex F, double theta, double phi);

        void check_weight();
        void check_transform();

    protected:
        fdcl_FFTS2_matrix_real nor_assoc_Legendre_poly(double cos_beta, int L);
        std::vector<double> compute_weight();
        double theta_k(int k);
        double phi_j(int j);

        fdcl_FFTS2_matrix_real nP;
        int L_4_check;

    private:
        fdcl_FFTS2_matrix_complex F_4_check;
        complex<double> f_4_check_transform(double theta, double phi);
};

class fdcl_FFTS2_real : public fdcl_FFTS2_complex
{
    public:
        fdcl_FFTS2_matrix_real y;
        fdcl_FFTS2_real(){};
        fdcl_FFTS2_real(int l_max);
        ~fdcl_FFTS2_real(){};
        void init(int l_max);

        fdcl_FFTS2_matrix_real spherical_harmonics(double theta, double phi, int L);
        fdcl_FFTS2_matrix_real forward_transform(std::function <double(double, double)>);
        double inverse_transform(fdcl_FFTS2_matrix_real F, double theta, double phi);

        fdcl_FFTSO3_matrix_complex T;
        fdcl_FFTSO3_matrix_complex matrix2rsph(int L);
        
        void check_transform();
    private:
        fdcl_FFTS2_matrix_real F_4_check;
        double f_4_check_transform(double theta, double phi);

};

#endif
