#include <iostream>
#include <vector>
#include <math.h> // pow
#include <iomanip>      // std::setprecision
#include <Eigen/Dense>


#include "fdcl_FFTSO3_matrix.hpp"
#include "fdcl_FFTSO3.hpp"
#include "fdcl_tictoc.hpp"
#include "misc_matrix_func.h"

using namespace std;
using namespace Eigen;

typedef Eigen::Matrix<double, 3, 3> Matrix3;

complex<double> myf(double a, double b, double g)
{
    return 1.0;
}

int main()
{
    int l_max=2;
    fdcl_FFTSO3_matrix_real d(l_max), d1(l_max), F_real(l_max);
    fdcl_FFTSO3_matrix_complex D(l_max), F0(l_max), F1(l_max), F(l_max);
    fdcl_FFTSO3 FFTSO3(l_max);
    fdcl_tictoc tictoc;
    double a, b, g;
    a=.12345;
    b=-0.234235;
    g=0.4324235;

   
//  FFTSO3.check_weight();
//  FFTSO3.check_wigner_d();

    // tictoc.tic();
    // F0=FFTSO3.forward_transform_0();
    // tictoc.toc("complex forward transform 0");

    // tictoc.tic();
    // F1=FFTSO3.forward_transform_1();
    // tictoc.toc("complex forward transform 1");
    
   // 
    // cout << FFTSO3.f(a,b,g) << endl;
    // tictoc.tic();
    // cout << FFTSO3.inverse_transform(F0,a,b,g) << endl;
    // tictoc.toc("complex inverse transform 0");
    // tictoc.tic();
    // cout << FFTSO3.inverse_transform(F1,a,b,g) << endl;
    // tictoc.toc("complex inverse transform 1");
    
    Vector3 eta1, eta2;
    Matrix3 R1, R2;
    eta1.setRandom();
    eta2.setRandom();
    
    R1=expm_SO3(eta1);
    R2=expm_SO3(eta2);

    FFTSO3.check_wigner_D_real();

    FFTSO3.check_forward_transform();

    // tictoc.tic();
    // F_real=FFTSO3.forward_transform_real();
    // tictoc.toc("real forward transform");
// 
    // tictoc.tic();
    // F_real=FFTSO3.forward_transform_real_0();
    // tictoc.toc("real forward transform 0");
// 
    // cout << FFTSO3.f_real(a,b,g) << endl;
// 
    // tictoc.tic();
    // cout << FFTSO3.inverse_transform(F_real,a,b,g) << endl; 
    // tictoc.toc("real inverse transform");
}
