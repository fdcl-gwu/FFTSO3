#include "fdcl_FFTS2.hpp"

fdcl::FFTS2_complex::FFTS2_complex(int l_max)
{
    init(l_max);
}

void fdcl::FFTS2_complex::init(int l_max)
{
    this->l_max=l_max;
    this->B=l_max+1;
    Y.init(l_max);
    F.init(l_max);
    weight.resize(2*B);
}

fdcl::FFTS2_matrix_complex fdcl::FFTS2_complex::spherical_harmonics(double theta, double phi, int L)
{
    fdcl::FFTS2_matrix_real nP(L);
    fdcl::FFTS2_matrix_complex Y(L);
    nP=nor_assoc_Legendre_poly(cos(theta),L);
    
#pragma omp parallel
{
    fdcl::omp_thread thr(omp_get_thread_num(),omp_get_num_threads());
    thr.range_closed(0,std::floor(0.5*L));

    for(int l=thr.i_init;l<=thr.i_term; l++)
    {
        for(int m=0; m<=l; m++)
            Y(l,m)=nP(l,m)*exp(I*(double)m*phi);

        for(int m=1; m<=l; m++)
            Y(l,-m)=pow(-1.,m)*std::conj(Y(l,m));

        int l_end=L-l;
        for(int m=0; m<=l_end; m++)
            Y(l_end,m)=nP(l_end,m)*exp(I*(double)m*phi);

        for(int m=1; m<=l_end; m++)
            Y(l_end,-m)=pow(-1.,m)*std::conj(Y(l_end,m));
    }
}

    return Y;
}

fdcl::FFTS2_matrix_real fdcl::FFTS2_complex::nor_assoc_Legendre_poly(double x, int L)
{
    // normalized associated Legendre polynomial of x with |x| < 1
    // Press, Teukolsky, Vetterling, Flannery, "Numerical Recepies", 3rd Edition, Sec 6.7 Spherical Harmonics
    fdcl::FFTS2_matrix_real nP(L);

    assert(abs(x) <= 1.0);

    nP.init(L);
    nP(0,0)=sqrt(1./(4.0*M_PI));

    double fact, oldfact, pmm, omx2;

    for(int l=1; l<=L; l++)
    {
        pmm=1.0; 

        omx2=(1.0-x)*(1.0+x);
        fact=1.0;
        for (int i=1;i<=l;i++) 
        {
            pmm *= omx2*fact/(fact+1.0);
            fact += 2.0;
        }
        pmm=sqrt((2*l+1)*pmm/(4.0*M_PI));

        if (l & 1) // change sign if abs(m) is odd: equivalent multiplying (-1)^m
            pmm=-pmm;

        nP(l,l)=pmm; // eqn (6.7.10)
    }

    for(int l=0; l<L; l++)
        nP(l+1,l)=x*sqrt(2.0*l+3.0)*nP(l,l);

    //#pragma omp for private (oldfact, fact, l, m)
    for(int m=0; m<=L; m++)
    {
        oldfact=sqrt(2.0*m+3.0);
        for(int l=m+2; l<=L; l++)
        { 
            fact=sqrt((4.0*l*l-1.0)/(l*l-m*m));
            nP(l,m)=(x*nP(l-1,m)-nP(l-2,m)/oldfact)*fact;
            oldfact=fact;
        }
    }

   for(int l=1; l<=L; l++)
        for(int m=1; m<=l; m++)
            nP(l,-m)=pow(-1.,m)*nP(l,m);
    
    return nP;
}

std::vector<double> fdcl::FFTS2_complex::compute_weight()
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

fdcl::FFTS2_matrix_complex fdcl::FFTS2_complex::forward_transform(std::function <complex<double>(double, double)> func)
{
    return forward_transform(func,0);
}

fdcl::FFTS2_matrix_complex fdcl::FFTS2_complex::forward_transform(std::function <complex<double>(double, double)> func, bool is_real)
{
    fdcl::FFTS2_matrix_complex F(l_max);
    F.setZero();
    compute_weight();

#pragma omp parallel
    {
        Eigen::VectorXcd F_k(2*B), func_k(2*B);
        Eigen::FFT<double> fft;
        double theta;

        fdcl::FFTS2_matrix_real nP_local(l_max);
        fdcl::FFTS2_matrix_complex F_local(l_max);

#pragma omp for private(theta)
        for(int k=0;k<2*B;k++)
        {
            theta=theta_k(k);
            for(int j=0; j<2*B; j++)
            {
                func_k(j)=func(theta,phi_j(j));
            }
            fft.fwd(F_k,func_k);
// 
            F_local.setZero();
            nP_local=nor_assoc_Legendre_poly(cos(theta),l_max);
            if(!is_real) // complex_valued func
            {
                for(int l=0; l<=l_max; l++)
                {    
                    for(int m=0; m<=l; m++)
                        F_local(l,m)+=weight[k]*nP_local(l,m)*F_k(m);
                    // 
                    for(int m=-l;  m<0; m++)
                        F_local(l,m)+=weight[k]*nP_local(l,m)*F_k(2*l_max+2+m);
                }
            }
            else // real-valued func
            {
                for(int l=0; l<=l_max; l++)
                    for(int m=0; m<=l; m++)
                        F_local(l,m)+=weight[k]*nP_local(l,m)*F_k(m);
            }
            // 

#pragma omp critical
            F=F+F_local;

        }
#pragma omp barrier
        
        if(is_real)
        {
#pragma omp for
            for(int l=0; l<=l_max; l++)
                for(int m=-l;  m<0; m++)
                    F(l,m)=pow(-1.,m)*std::conj(F(l,-m));
        }

    }
    return F;
}

complex<double> fdcl::FFTS2_complex::inverse_transform(fdcl::FFTS2_matrix_complex F, double theta, double phi)
{
    complex<double> y={0., 0.};
    fdcl::FFTS2_matrix_complex Y(F.l_max);
    Y=spherical_harmonics(theta,phi,F.l_max);

// required OpenMp > 4.0
// #pragma omp declare reduction \
    // (complex_sum : std::complex<double> : omp_out+=omp_in) \
    // initializer(omp_priv={0.,0.})
// 
// #pragma omp parallel for reduction(complex_sum:y)
    // for(int l=0; l<=l_max; l++)
        // for(int m=-l; m<=l; m++)
            // y+=F(l,m)*Y(l,m);

#pragma omp parallel
    {
        fdcl::omp_thread thr(omp_get_thread_num(),omp_get_num_threads());
        thr.range_closed(0,F.l_max);
        std::complex<double> y_local={0.,0.};

        for(int l=thr.i_init; l<=thr.i_term; l++)
        {
            for(int m=0; m<=l; m++)
                y_local+=F(l,m)*Y(l,m);
            int l_end = F.l_max-l;
            for(int m=-l_end; m<0; m++)
                y_local+=F(l_end,m)*Y(l_end,m);
        }
#pragma omp critical
        y+=y_local;

    }

    return y;
}

complex<double> fdcl::FFTS2_complex::inverse_transform(double theta, double phi)
{
    return inverse_transform(this->F,theta,phi);
}

double fdcl::FFTS2_complex::check_weight()
{
    fdcl::FFTS2_matrix_complex Y(2*B-1);
    std::vector<complex<double>> sum;
    double error=1.0;
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
    
    if(check_verbose)
    {
        cout << "fdcl::FFTS2_complex::check_weight" << endl;
        cout << "\\sum_k w_k Y^l_m(\\theta_k,0) * B / \\sqrt{\\pi} = \\delta_{0,l}" << endl; 
        for (int l=0;l<2*B;l++)
            cout << "l=" << l << ": " << sum[l]*(double)B/sqrt(M_PI) << endl;
    }
    error = abs(sum[0]*(double)B/sqrt(M_PI)-1.0);
    for (int l=1;l<2*B;l++)
        error+= abs(sum[l]*(double)B/sqrt(M_PI));

    int j, k, l;
    fdcl::FFTS2_matrix_complex Delta(2*B-1);
    Delta.setZero();
    for(k=0;k<2*B;k++)
        for(j=0;j<2*B;j++)
            for(l=0;l<2*B;l++)
                Delta[l]+=weight[k]*spherical_harmonics(theta_k(k),phi_j(j),2*B-1)[l];
    
    if(check_verbose)
    {
        cout << "\\sum_{j,k} w_k Y(theta_k, phi_j) / (2\\sqrt{\\pi}) = \\delta_{l,0}\\delta_{m,0}" << endl;
        for (int l=0;l<2*B;l++)
            cout << "l=" << l << ": " << Delta[l].norm()/(2.*sqrt(M_PI)) << endl;
    }
    error += abs(Delta(0,0)/(2.*sqrt(M_PI))-1.0);
    for (int l=1;l<2*B;l++)
        error+=abs(Delta[l].norm()/(2.*sqrt(M_PI)));

    cout << "fdcl::FFTS2_complex::check_weight: error = " << error << endl;
    return error;
}

void fdcl::FFTS2_complex::check_all()
{
    check_weight();
    check_transform();
    cout << endl;
}

double fdcl::FFTS2_complex::theta_k(int k)
{
    return ((double)(2*k+1))*M_PI/4./((double)B);
}

double fdcl::FFTS2_complex::phi_j(int j)
{
    return ((double)j)*M_PI/((double)B);
}

double fdcl::FFTS2_complex::check_transform()
{
    double error=1.0;

    F_4_check.init(l_max);
    F_4_check.setRandom();

    auto func= std::bind(&fdcl::FFTS2_complex::f_4_check_transform, this, std::placeholders::_1, std::placeholders::_2);

    error = (F_4_check-forward_transform(func)).norm();

    cout << "fdcl::FFTS2_complex::check_transform: l_max=" << l_max << " : error = " << error << endl;

    return error;
}

complex<double> fdcl::FFTS2_complex::f_4_check_transform(double theta, double phi)
{
    return inverse_transform(F_4_check, theta, phi);
}


fdcl::FFTSO3_matrix_complex fdcl::FFTS2_real::matrix2rsph(int L)
{
    double one_over_sqrt_2=1./sqrt(2.);
    T.init(L);

#pragma omp parallel
    {
        fdcl::omp_thread thr(omp_get_thread_num(),omp_get_num_threads());
        thr.range_closed(0,L);
        int l,m;

        for(l=thr.i_init;l<=thr.i_term;l++)
        {
            T(l,0,0)=1.;
            for(m=1;m<=l;m++)
            {
                T(l,m,-m)=one_over_sqrt_2;
                T(l,m,m)=pow(-1.,m)*one_over_sqrt_2;
            }
            int l_end=l_max-l;
            for(m=-l_end;m<0;m++)
            {
                T(l_end,m,m)=I*one_over_sqrt_2;
                T(l_end,m,-m)=-I*one_over_sqrt_2*pow(-1.,m);
            }
        }

    }

    // for(l=0;l<=L;l++)
    // {
       // for(m=-l;m<0;m++)
        // {
            // T(l,m,m)=I/sqrt(2.);
            // T(l,m,-m)=-I/sqrt(2.)*pow(-1.,m);
        // }
        // T(l,0,0)=1.;
        // for(m=1;m<=l;m++)
        // {
            // T(l,m,-m)=1./sqrt(2.);
            // T(l,m,m)=pow(-1.,m)/sqrt(2.);
        // }
    // }
// 
    return T;
}

fdcl::FFTS2_real::FFTS2_real(int l_max)
{
    init(l_max);
}

void fdcl::FFTS2_real::init(int l_max)
{
    fdcl::FFTS2_complex::init(l_max);
    y.init(l_max);
    F.init(l_max);
}

fdcl::FFTS2_matrix_real fdcl::FFTS2_real::spherical_harmonics(double theta, double phi, int L)
{
    fdcl::FFTS2_matrix_real nP(L), y(L);

    nP=nor_assoc_Legendre_poly(cos(theta),L);
#pragma omp parallel
    {
        double tmp;
        fdcl::omp_thread thr(omp_get_thread_num(),omp_get_num_threads());
        thr.range_closed(0,std::floor(0.5*L));

        for(int l=thr.i_init; l<=thr.i_term; l++)
        {
            y(l,0)=nP(l,0);
            for(int m=1; m<=l; m++)
            {
                tmp=sqrt(2.)*pow(-1.,m)*nP(l,m);
                y(l,m)=tmp*cos((double)m*phi);
                y(l,-m)=tmp*sin((double)m*phi);
            }
            int l_end=L-l;
            y(l_end,0)=nP(l_end,0);
            for(int m=1; m<=l_end; m++)
            {
                tmp=sqrt(2.)*pow(-1.,m)*nP(l_end,m);
                y(l_end,m)=tmp*cos((double)m*phi);
                y(l_end,-m)=tmp*sin((double)m*phi);
            }
        }
    }

    // alternative method: conversion from complex harmonics: slower
    // y.init(L);
    // fdcl::FFTS2_complex::spherical_harmonics(theta,phi,L);
    // matrix2rsph(L);
    // for(int l=0; l<=L; l++)
       // y[l]= (T[l]*Y[l]).real();

    return y;
}

fdcl::FFTS2_matrix_real fdcl::FFTS2_real::forward_transform(std::function <double(double, double)> func)
{
    fdcl::FFTS2_matrix_complex F_complex(l_max);

    fdcl::FFTSO3_matrix_complex T(l_max);
    F_complex=fdcl::FFTS2_complex::forward_transform(func,1);
    T=matrix2rsph(l_max);

#pragma omp parallel for
    for(int l=0; l<=l_max;l++)
        F[l]=(T[l].conjugate()*F_complex[l]).real();

    return F;

    // alternative method: forward transform with RSH: slower
    // fdcl::FFTS2_matrix_real F(l_max);
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


void fdcl::FFTS2_real::check_all()
{
    check_transform();
    cout << endl;
}

double fdcl::FFTS2_real::check_transform()
{
    double error=1.0;
    F_4_check.init(l_max);
    F_4_check.setRandom();

    auto func= std::bind(&fdcl::FFTS2_real::f_4_check_transform, this, std::placeholders::_1, std::placeholders::_2);

    error = (F_4_check-forward_transform(func)).norm();
    cout << "fdcl::FFTS2_real::check_transform: l_max=" << l_max  << " : error = " << error << endl;

    return error;
}

double fdcl::FFTS2_real::f_4_check_transform(double theta, double phi)
{
    return inverse_transform(F_4_check, theta, phi);
}

double fdcl::FFTS2_real::inverse_transform(fdcl::FFTS2_matrix_real F, double theta, double phi)
{
    double z=0.;
    fdcl::FFTS2_matrix_real y(F.l_max);
    y=spherical_harmonics(theta,phi,F.l_max);
    
#pragma omp parallel for reduction(+:z)
    for(int l=0; l<=F.l_max; l++)
        for(int m=-l; m<=l; m++)
            z+=F(l,m)*y(l,m);

    return z;
}

double fdcl::FFTS2_real::inverse_transform(double theta, double phi)
{
    return inverse_transform(this->F, theta, phi);
}

