#include "fdcl_FFTS2.hpp"

fdcl_FFTS2_complex::fdcl_FFTS2_complex(int l_max)
{
    init(l_max);
}

void fdcl_FFTS2_complex::init(int l_max)
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
    return forward_transform(func,0);
}

fdcl_FFTS2_matrix_complex fdcl_FFTS2_complex::forward_transform(std::function <complex<double>(double, double)> func, bool is_real)
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

    // alternative method without fft: much slower
    // F_km_2.resize(2*B,2*l_max+1);
    // F_km_2.setZero();
    // for(int k=0; k<2*B; k++)
    // {
        // for(int j=0; j<2*B; j++)
        // {
            // for(int m=-l_max; m<=l_max; m++)
            // {
                // F_km_2(k,m+l_max)+=func(theta_k(k),phi_j(j))*exp(-I*(double)m*phi_j(j));
            // }
        // }
    // }
 
    F.setZero();
    compute_weight();

    if(is_real)
    {
        for(int k=0;k<2*B;k++)
        {
            nor_assoc_Legendre_poly(cos(theta_k(k)),l_max);
            for(int l=0; l<=l_max; l++)
                for(int m=0; m<=l; m++)
                    F(l,m)+=weight[k]*nP(l,m)*F_km(k,m);
        }

        for(int l=0; l<=l_max; l++)
            for(int m=-l;  m<0; m++)
                F(l,m)=pow(-1.,m)*std::conj(F(l,-m));
    }
    else
    {
        for(int k=0;k<2*B;k++)
        {
            nor_assoc_Legendre_poly(cos(theta_k(k)),l_max);
            for(int l=0; l<=l_max; l++)
            {    
                for(int m=0; m<=l; m++)
                    F(l,m)+=weight[k]*nP(l,m)*F_km(k,m);

                for(int m=-l;  m<0; m++)
                    F(l,m)+=weight[k]*nP(l,m)*(F_km(k,2*l_max+2+m));
            }
        }
    }

    return F;
}

complex<double> fdcl_FFTS2_complex::inverse_transform(fdcl_FFTS2_matrix_complex F, double theta, double phi)
{
    complex<double> y={0., 0.};
    init(F.l_max);
    spherical_harmonics(theta,phi,F.l_max);
    
    for(int l=0; l<=F.l_max; l++)
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

void fdcl_FFTS2_complex::check_transform()
{
    F_4_check.init(l_max);
    F_4_check.setRandom();

    auto func= std::bind(&fdcl_FFTS2_complex::f_4_check_transform, this, std::placeholders::_1, std::placeholders::_2);

    cout << "fdcl_FFTS2_complex::check_transform: l_max=" << l_max << endl;
    cout << "error = " << (F_4_check-forward_transform(func)).norm() << endl;

    init(l_max);
}

complex<double> fdcl_FFTS2_complex::f_4_check_transform(double theta, double phi)
{
    return inverse_transform(F_4_check, theta, phi);
}


fdcl_FFTSO3_matrix_complex fdcl_FFTS2_real::matrix2rsph(int L)
{
    int l,m;
    T.init(L);

    for(l=0;l<=L;l++)
    {
       for(m=-l;m<0;m++)
        {
            T(l,m,m)=I/sqrt(2.);
            T(l,m,-m)=-I/sqrt(2.)*pow(-1.,m);
        }
        T(l,0,0)=1.;
        for(m=1;m<=l;m++)
        {
            T(l,m,-m)=1./sqrt(2.);
            T(l,m,m)=pow(-1.,m)/sqrt(2.);
        }
    }

    return T;
}

fdcl_FFTS2_real::fdcl_FFTS2_real(int l_max)
{
    init(l_max);
}

void fdcl_FFTS2_real::init(int l_max)
{
    fdcl_FFTS2_complex::init(l_max);
    y.init(l_max);
}

fdcl_FFTS2_matrix_real fdcl_FFTS2_real::spherical_harmonics(double theta, double phi, int L)
{
    y.init(L);
    double tmp;

    nor_assoc_Legendre_poly(cos(theta),L);
    for(int l=0; l<=L; l++)
    {
        y(l,0)=nP(l,0);
        for(int m=1; m<=l; m++)
        {
            tmp=sqrt(2.)*pow(-1.,m)*nP(l,m);
            y(l,m)=tmp*cos((double)m*phi);
            y(l,-m)=tmp*sin((double)m*phi);
        }
    }

    // alternative method: conversion from complex harmonics: slower
    // y.init(L);
    // fdcl_FFTS2_complex::spherical_harmonics(theta,phi,L);
    // matrix2rsph(L);
    // for(int l=0; l<=L; l++)
       // y[l]= (T[l]*Y[l]).real();

    return y;
}

fdcl_FFTS2_matrix_real fdcl_FFTS2_real::forward_transform(std::function <double(double, double)> func)
{
    fdcl_FFTS2_matrix_complex F(l_max);
    fdcl_FFTS2_matrix_real F_out(l_max);

    fdcl_FFTSO3_matrix_complex T(l_max);
    F=fdcl_FFTS2_complex::forward_transform(func,1);
    T=matrix2rsph(l_max);
    for(int l=0; l<=l_max;l++)
        F_out[l]=(T[l].conjugate()*F[l]).real();

    return F_out;

    // alternative method: forward transform with RSH: slower
    // fdcl_FFTS2_matrix_real F(l_max);
    // Eigen::Matrix<double, Dynamic, Dynamic> F_km;
    // double tmp;
    // int m;
// 
    // F_km.resize(2*B,2*l_max+1);
    // F_km.setZero();
    // for(int k=0; k<2*B; k++)
        // for(int j=0; j<2*B; j++)
        // {
            // tmp=func(theta_k(k),phi_j(j));
            // m=0;
            // F_km(k,m+l_max)+=tmp;
// 
            // for(m=1; m<=l_max; m++)
            // {
                // F_km(k,m+l_max)+=tmp*sqrt(2.)*pow(-1.,m)*cos((double)m*phi_j(j));
                // F_km(k,-m+l_max)+=tmp*sqrt(2.)*pow(-1.,m)*sin((double)m*phi_j(j));
            // }
        // }
// 
    // F.setZero();
    // compute_weight();
// 
    // for(int k=0;k<2*B;k++)
    // {
        // nor_assoc_Legendre_poly(cos(theta_k(k)),l_max);
        // for(int l=0; l<=l_max; l++)
            // for(int m=-l; m<=l; m++)
                // F(l,m)+=weight[k]*nP(l,abs(m))*F_km(k,m+l_max);
    // }
// 
    // return F;
}

void fdcl_FFTS2_real::check_transform()
{
    F_4_check.init(l_max);
    F_4_check.setRandom();

    auto func= std::bind(&fdcl_FFTS2_real::f_4_check_transform, this, std::placeholders::_1, std::placeholders::_2);

    cout << "fdcl_FFTS2_real::check_transform: l_max=" << l_max << endl;
    cout << "error = " << (F_4_check-forward_transform(func)).norm() << endl;

    init(l_max);
}

double fdcl_FFTS2_real::f_4_check_transform(double theta, double phi)
{
    return inverse_transform(F_4_check, theta, phi);
}

double fdcl_FFTS2_real::inverse_transform(fdcl_FFTS2_matrix_real F, double theta, double phi)
{
    double z=0.;
    init(F.l_max);
    spherical_harmonics(theta,phi,F.l_max);
    
    for(int l=0; l<=F.l_max; l++)
        for(int m=-l; m<=l; m++)
            z+=F(l,m)*y(l,m);

    return z;
}

