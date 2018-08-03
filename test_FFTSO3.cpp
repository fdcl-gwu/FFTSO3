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

class fdcl_FFTS2_complex
{
    public:
       int l_max;
       fdcl_FFTS2_matrix_complex Y;

       fdcl_FFTS2_complex(){};
       fdcl_FFTS2_complex(int l_max);
       ~fdcl_FFTS2_complex(){};

       fdcl_FFTS2_matrix_complex spherical_harmonics(double theta, double phi, int L);

       fdcl_FFTS2_matrix_real nor_assoc_Legendre_poly(double cos_beta, int L);
    private:
       fdcl_FFTS2_matrix_real N, nP;
};

fdcl_FFTS2_complex::fdcl_FFTS2_complex(int l_max)
{
    this->l_max=l_max;
    Y.init(l_max);
}

fdcl_FFTS2_matrix_complex fdcl_FFTS2_complex::spherical_harmonics(double theta, double phi, int L)
{
    nor_assoc_Legendre_poly(cos(theta),L);
    Y.init(L);

    for(int l=0;l<=L; l++)
    {
        for(int m=0; m<=l; m++)
            Y(l,m)=nP(l,m)*exp(I*(double)m*phi);

        for(int m=1; m<=l; m++)
            Y(l,-m)=pow(-1.,m)*std::conj(Y(l,m));
    }

    return Y;
}


fdcl_FFTS2_matrix_real fdcl_FFTS2_complex::nor_assoc_Legendre_poly(double x, int L)
{
    // normalized associated Legendre polynomial of x with |x| < 1
    // Press, Teukolsky, Vetterling, Flannery, "Numerical Recepies", 3rd Edition, Sec 6.7 Spherical Harmonics
    int i, m, l;
    double fact, oldfact, pmm, omx2;

    nP.init(L);

    assert(abs(x) <= 1.0);

    nP(0,0)=sqrt(1./(4.0*M_PI));

    for(l=1; l<=L; l++)
    {
        pmm=1.0; 

        omx2=(1.0-x)*(1.0+x);
        fact=1.0;
        for (i=1;i<=l;i++) 
        {
            pmm *= omx2*fact/(fact+1.0);
            fact += 2.0;
        }
        pmm=sqrt((2*l+1)*pmm/(4.0*M_PI));

        if (l & 1) // change sign if abs(m) is odd: equivalent multiplying (-1)^m
            pmm=-pmm;

        nP(l,l)=pmm; // eqn (6.7.10)
    }

    for(l=0; l<L; l++)
        nP(l+1,l)=x*sqrt(2.0*l+3.0)*nP(l,l);

    for(m=0; m<=L; m++)
    {
        oldfact=sqrt(2.0*m+3.0);
        for(l=m+2; l<=L; l++)
        { 
            fact=sqrt((4.0*l*l-1.0)/(l*l-m*m));
            nP(l,m)=(x*nP(l-1,m)-nP(l-2,m)/oldfact)*fact;
            oldfact=fact;
        }
    }

    return nP;
}



int main()
{
    int l_max=3;
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
    cout << FFTS2.spherical_harmonics(0.12345,9.87654,3) << endl;

}
