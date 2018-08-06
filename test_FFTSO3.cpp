#include <iostream>
#include <vector>
#include <math.h> // pow
#include <iomanip>      // std::setprecision
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

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
       int l_max, B;
       fdcl_FFTS2_matrix_complex Y;
       std::vector<double> weight;

       fdcl_FFTS2_complex(){};
       fdcl_FFTS2_complex(int l_max);
       ~fdcl_FFTS2_complex(){};

       fdcl_FFTS2_matrix_complex spherical_harmonics(double theta, double phi, int L);


       fdcl_FFTS2_matrix_complex forward_transform(std::function <complex<double>(double, double)>);
       complex<double> inverse_transform(fdcl_FFTS2_matrix_complex F, double theta, double phi);

       void check_weight();

    private:
       fdcl_FFTS2_matrix_real nor_assoc_Legendre_poly(double cos_beta, int L);
       std::vector<double> compute_weight();
       double theta_k(int k);
       double phi_j(int j);
    private:
       fdcl_FFTS2_matrix_real nP;
};

fdcl_FFTS2_complex::fdcl_FFTS2_complex(int l_max)
{
    this->l_max=l_max;
    this->B=l_max+1;
    Y.init(l_max);
    weight.resize(2*B);
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

    // copy for strictly negative m
    for(l=1; l<=L; l++)
        for(m=1; m<=l; m++)
            nP(l,-m)=pow(-1.,m)*nP(l,m);

    return nP;
}


std::vector<double> fdcl_FFTS2_complex::compute_weight()
{   
    int j, k;
    double factor;
    double sum;

    factor = M_PI/((double)(4*B)) ;

    for(j=0;j<2*B;j++)
    {
        sum=0.0;
        for(k=0;k<B;k++)
            sum+=1./((double)(2*k+1))*sin((double)((2*j+1)*(2*k+1))*factor);
        
        sum*=2.*M_PI/((double)B*B)*sin((double)(2*j+1)*factor);
      
        weight[j]=sum;
    }

    return weight;
}

fdcl_FFTS2_matrix_complex fdcl_FFTS2_complex::forward_transform(std::function <complex<double>(double, double)> func)
{
    fdcl_FFTS2_matrix_complex F(l_max);
    Eigen::Matrix<complex<double>, Dynamic, Dynamic> F_km, F_km_2;
    Eigen::VectorXcd func_k(2*B), tmp_out(2*B);
    Eigen::FFT<double> fft;
    double theta;

    F_km.resize(2*B,2*B);

    for(int k=0;k<2*B;k++)
    {
        theta=theta_k(k);
        for(int j=0; j<2*B; j++)
            func_k(j)=func(theta,phi_j(j));
        fft.fwd(tmp_out,func_k);

        F_km.row(k)=tmp_out;
    }

    F_km_2.resize(2*B,2*l_max+1);
    F_km_2.setZero();
    for(int k=0; k<2*B; k++)
    {
        for(int j=0; j<2*B; j++)
        {
            for(int m=-l_max; m<=l_max; m++)
            {
                F_km_2(k,m+l_max)+=func(theta_k(k),phi_j(j))*exp(-I*(double)m*phi_j(j));
            }
        }
    }

    F.setZero();
    compute_weight();
    for(int k=0;k<2*B;k++)
    {
        nor_assoc_Legendre_poly(cos(theta_k(k)),l_max);
        for(int l=0; l<=l_max; l++)
        {    
            for(int m=-l; m<=l; m++)
                F(l,m)+=weight[k]*nP(l,m)*F_km_2(k,l_max+m);
        }
    }

    return F;
}

complex<double> fdcl_FFTS2_complex::inverse_transform(fdcl_FFTS2_matrix_complex F, double theta, double phi)
{
    complex<double> y={0., 0.};
    spherical_harmonics(theta,phi,l_max);
    
    for(int l=0; l<=l_max; l++)
        for(int m=-l; m<=l; m++)
            y+=F(l,m)*Y(l,m);

    return y;
}


void fdcl_FFTS2_complex::check_weight()
{
    fdcl_FFTS2_matrix_complex Y(2*B-1);
    std::vector<complex<double>> sum;
    sum.resize(2*B);
    for (int l=0;l<2*B;l++)
        sum[l]=0.;

    this->compute_weight();
    
    for (int k=0;k<2*B; k++)
    {
        Y=spherical_harmonics(theta_k(k),0.,2*B-1);
        for(int l=0;l<2*B;l++)
        {
            sum[l]+=Y(l,0)*weight[k];
        }
    }
    
    cout << "fdcl_FFTS2_complex::check_weight" << endl;
    cout << "\\sum_k w_k Y^l_m(\\theta_k,0) * B / \\sqrt{\\pi} = \\delta_{0,l}" << endl; 
    for (int l=0;l<2*B;l++)
        cout << "l=" << l << ": " << sum[l]*(double)B/sqrt(M_PI) << endl;

    int j, k, l;
    fdcl_FFTS2_matrix_complex Delta(2*B-1);
    Delta.setZero();
    for(k=0;k<2*B;k++)
        for(j=0;j<2*B;j++)
            for(l=0;l<2*B;l++)
                Delta[l]+=weight[k]*spherical_harmonics(theta_k(k),phi_j(j),2*B-1)[l];
    
    cout << "\\sum_{j,k} w_k Y(theta_k, phi_j) / (2\\sqrt{\\pi}) = \\delta_{l,0}\\delta_{m,0}" << endl;
    for (int l=0;l<2*B;l++)
        cout << "l=" << l << ": " << Delta[l].norm()/(2.*sqrt(M_PI)) << endl;

    
}

double fdcl_FFTS2_complex::theta_k(int k)
{
    return ((double)(2*k+1))*M_PI/4./((double)B);
}

double fdcl_FFTS2_complex::phi_j(int j)
{
    return ((double)j)*M_PI/((double)B);
}

complex<double> myf_S2(double theta, double phi)
{
    fdcl_FFTS2_complex FFTS2;
    FFTS2.spherical_harmonics(theta,phi,3);

    return cos(theta)*cos(phi)+cos(theta)*sin(phi);
    // return FFTS2.Y(1,1);
    // return 1.0;
}

int main()
{
    int l_max=80;
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
    fdcl_FFTS2_matrix_complex F_SH;
    F_SH=FFTS2.forward_transform(myf_S2);

    double theta=0.12345, phi=0.8765432;
    cout << FFTS2.inverse_transform(F_SH, theta, phi) << endl;
    cout << myf_S2(theta, phi) << endl;
// 

    // Eigen::FFT<double> fft;
// 
    // std::vector<double> in;
    // in.resize(10);
    // for(int i=0; i<10; i++)
        // in[i]=sin((double)i);
// 
    // std::vector<complex<double>> out;
    // fft.fwd(out, in);
    // 
    // cout << out.size() << endl;
    // for (int i=0; i<10; i++)
        // cout << out[i] << endl;
// 
}
