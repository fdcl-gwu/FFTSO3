#include <iostream>
#include <vector>
#include <math.h> // pow
#include <iomanip>      // std::setprecision
#include <Eigen/Dense>
#include "fdcl_FFTSO3_matrix.hpp"
#include "fdcl_FFTSO3.hpp"
#include "misc_matrix_func.h"

using namespace std;
using namespace Eigen;

typedef Eigen::Matrix<double, 3, 3> Matrix3;

class fdcl_tictoc
{
public:
    std::chrono::steady_clock::time_point t0, t1;
    fdcl_tictoc(){};
    void tic();
    void toc();
};

void fdcl_tictoc::tic()
{
    t0 = std::chrono::steady_clock::now();  
}
void fdcl_tictoc::toc()
{
    t1 = std::chrono::steady_clock::now();  
    cout << "fdcl_tictoc = " << std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()/1e6 << " sec" << endl;
}

int main()
{
    int l_max=2;
    fdcl_FFTSO3_matrix_real d(l_max), d1(l_max);
    fdcl_FFTSO3_matrix_complex D(l_max), F0(l_max), F1(l_max);
    fdcl_FFTSO3 FFTSO3(l_max);
    fdcl_tictoc tictoc;
    
//  FFTSO3.check_weight();
//  FFTSO3.check_wigner_d();

    // tictoc.tic();
    // F0=FFTSO3.forward_transform_0();
    // tictoc.toc();
// 
    // cout << F0 << endl;
    // tictoc.tic();
    // F1=FFTSO3.forward_transform_1();
    // tictoc.toc();
// 
    // 
    // for(int l=0;l<=l_max;l++)
        // cout << (F0[l]-F1[l]).norm() << endl;
    // 
    // cout << FFTSO3.f(1.,2.,3.) << endl;
    // cout << FFTSO3.inverse_transform(F0,1.,2.,3.) << endl;
    // cout << FFTSO3.inverse_transform(F1,1.,2.,3.) << endl;
    // 
    Vector3 eta1, eta2;
    Matrix3 R1, R2;
    eta1.setRandom();
    eta2.setRandom();
    
    R1=expm_SO3(eta1);
    R2=expm_SO3(eta2);

    double a, b, g;
    a=.12345;
    b=-9.234235;
    g=0.4324235;

    tictoc.tic();
    for(int i=0;i<100;i++)
        FFTSO3.wigner_D_real(a,b,g,l_max);
    tictoc.toc();

    tictoc.tic();
    for(int i=0;i<100;i++)
        FFTSO3.wigner_D_real_converted(a,b,g,l_max);
    tictoc.toc();

    tictoc.tic();
    for(int i=0;i<100;i++)
        FFTSO3.wigner_D_real_Phi(a,b,g,l_max);
    tictoc.toc();

}
