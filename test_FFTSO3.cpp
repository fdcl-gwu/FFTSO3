#include <iostream>
#include <vector>
#include <math.h> // pow
#include <iomanip>      // std::setprecision
#include <Eigen/Dense>

#include "fdcl_FFTSO3_matrix.hpp"
#include "fdcl_FFTSO3_complex.hpp"
#include "fdcl_FFTSO3_real.hpp"
#include "fdcl_tictoc.hpp"
#include "misc_matrix_func.h"

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

int main()
{
    int l_max=2;
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
    FFTSO3.check_Clebsch_Gordon();
    

    // FFTSO3_real.check_weight();
    // FFTSO3_real.check_wigner_d();
    // FFTSO3_real.check_wigner_D_real();
    // FFTSO3_real.check_forward_transform();
    FFTSO3_real.check_Clebsch_Gordon();
}
