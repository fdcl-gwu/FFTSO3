#include <iostream>
#include <vector>
#include <math.h> // pow
#include <iomanip>      // std::setprecision
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

#include "fdcl_FFTSO3_matrix.hpp"
#include "fdcl_FFTSO3_complex.hpp"
#include "fdcl_FFTSO3_real.hpp"
#include "fdcl_FFTS2.hpp"
#include "fdcl_tictoc.hpp"
#include "misc_matrix_func.h"
#include "fdcl_FFTS2_matrix.hpp"

using namespace std;
using namespace Eigen;

typedef Eigen::Matrix<double, 3, 3> Matrix3;

double myrf(double a, double b, double g)
{
    return Euler3232R(a,b,g).trace();
}

double myrfR(Matrix3 R)
{
    return R.trace();
}

complex<double> myf(double a, double b, double g)
{
    return Euler3232R(a,b,g).trace();
}

complex<double> myfR(Matrix3 R)
{
    return R.trace();
}

complex<double> myf_S2(double theta, double phi)
{
    fdcl_FFTS2_complex FFTS2;
    FFTS2.spherical_harmonics(theta,phi,3);

    return cos(theta)*cos(phi)+cos(theta)*sin(phi);
    // return FFTS2.Y(1,1);
    // return 1.0;
}
double myrf_S2(double theta, double phi)
{
    fdcl_FFTS2_complex FFTS2;
    FFTS2.spherical_harmonics(theta,phi,3);

    return cos(theta)*cos(phi)+cos(theta)*sin(phi);
    // return FFTS2.Y(1,1);
    // return 1.0;
}

int main()
{
    int l_max=5;
    fdcl_FFTSO3_matrix_real d(l_max), d1(l_max), F_real(l_max);
    fdcl_FFTSO3_matrix_complex D(l_max), F0(l_max), F1(l_max), F(l_max);
    fdcl_FFTSO3_complex FFTSO3(l_max);
    fdcl_FFTSO3_real FFTSO3_real(l_max);
    fdcl_tictoc tictoc;
    double a, b, g;
    a=.12345;
    b=-0.234235;
    g=0.4324235;
    Matrix3 R;
    R.setIdentity();
    Vector3 eta1, eta2;
    Matrix3 R1, R2;
    eta1.setRandom();
    eta2.setRandom();
    
    R1=expm_SO3(eta1);
    R2=expm_SO3(eta2);
    
    // FFTSO3.check_weight();
    // FFTSO3.check_wigner_d();
    // FFTSO3.check_forward_transform();
    // FFTSO3.check_Clebsch_Gordon();
    
    // FFTSO3_real.check_weight();
    // FFTSO3_real.check_wigner_d();
    // FFTSO3_real.check_wigner_D_real();
    // FFTSO3_real.check_forward_transform();
    // FFTSO3_real.check_Clebsch_Gordon();
    // FFTSO3_real.check_deriv_U();
    //

    fdcl_FFTS2_complex FFTS2(l_max);
    // FFTS2.check_weight();
    FFTS2.check_transform();

    fdcl_FFTS2_real RFFTS2(l_max);
    RFFTS2.check_transform();

}
