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
    L_4_check=100;
    F_4_check.init(L_4_check);
    F_4_check.setRandom();
    init(L_4_check);

    auto func= std::bind(&fdcl_FFTS2_complex::f_4_check_transform, this, std::placeholders::_1, std::placeholders::_2);

    cout << "fdcl_FFTS2_complex::check_transform: l_max=" << L_4_check << endl;
    cout << "error = " << (F_4_check-forward_transform(func)).norm() << endl;
}

complex<double> fdcl_FFTS2_complex::f_4_check_transform(double theta, double phi)
{
    return inverse_transform(F_4_check, theta, phi);
}


