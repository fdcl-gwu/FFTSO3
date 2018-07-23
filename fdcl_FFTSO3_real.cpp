#include "fdcl_FFTSO3_matrix.hpp"
#include "fdcl_FFTSO3_real.hpp"

fdcl_FFTSO3_real::fdcl_FFTSO3_real(int l_max)
{
    this->l_max=l_max;  
    this->B=l_max+1;
    d_beta.resize(2*B);
    weight.resize(2*B);
}

Eigen::VectorXd fdcl_FFTSO3_real::Legendre_poly(double x, int N)
{
    // return the value of Legendre Polynomial of x upto the order N
    Eigen::VectorXd P;
    int n;
    
    P.resize(N+1);
    
    P(0)=1.;
    if (N>=1)
        P(1)=x;
    
    for(n=1;n<N;n++)
        P(n+1)=1./((double)n+1)*( ((double)2*n+1)*x*P(n)-((double)n)*P(n-1) );
    
    return P;
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_real::wigner_d_explicit(double beta)
{
    // D. Varshalovich, A. Moskalev, and V. Khersonskii, Quantum Theory of Angular Momentum, World Scientific, 1988, Chapter 4
    
    fdcl_FFTSO3_matrix<double> d(3);
    double cb, sb, sb2, cb2, tb2;
    cb=cos(beta);
    sb=sin(beta);
    sb2=sin(beta/2.);
    cb2=cos(beta/2.);
    tb2=tan(beta/2.);
    
    d(0,0,0)=1.;
    
    d(1,1,1)=(1.+cb)/2.;
    d(1,1,0)=-sb/sqrt(2.);
    d(1,1,-1)=(1.-cb)/2.;
    d(1,0,1)=sb/sqrt(2.);
    d(1,0,0)=cb;
    d(1,0,-1)=-sb/sqrt(2.);
    d(1,-1,1)=(1-cb)/2.;
    d(1,-1,0)=sb/sqrt(2.);
    d(1,-1,-1)=(1.+cb)/2.;
        
    d(2,2,2)=pow(1+cb,2)/4.;
    d(2,2,1)=-sb*(1.+cb)/2.;
    d(2,2,0)=1./2.*sqrt(3./2.)*pow(sb,2);
    d(2,2,-1)=-sb*(1.-cb)/2.;
    d(2,2,-2)=pow(1-cb,2)/4.;
    d(2,1,2)=sb*(1.+cb)/2.;
    d(2,1,1)=(2.*pow(cb,2)+cb-1.)/2.;
    d(2,1,0)=-sqrt(3./2.)*sb*cb;
    d(2,1,-1)=-(2.*pow(cb,2)-cb-1.)/2.;
    d(2,1,-2)=-sb*(1-cb)/2.;
    d(2,0,2)=1./2.*sqrt(3./2.)*pow(sb,2);
    d(2,0,1)=sqrt(3./2.)*sb*cb;
    d(2,0,0)=(3.*pow(cb,2)-1.)/2.;
    d(2,0,-1)=-sqrt(3./2.)*sb*cb;
    d(2,0,-2)=1./2.*sqrt(3./2.)*pow(sb,2);
    d(2,-1,2)=sb*(1.-cb)/2.;
    d(2,-1,1)=-(2.*pow(cb,2)-cb-1.)/2.;
    d(2,-1,0)=sqrt(3./2.)*sb*cb;
    d(2,-1,-1)=(2*pow(cb,2)+cb-1.)/2.;
    d(2,-1,-2)=-sb*(1+cb)/2.;
    d(2,-2,2)=pow(1.-cb,2)/4.;
    d(2,-2,1)=sb*(1.-cb)/2.;
    d(2,-2,0)=1./2.*sqrt(3./2.)*pow(sb,2);
    d(2,-2,-1)=sb*(1.+cb)/2.;
    d(2,-2,-2)=pow(1+cb,2)/4.;
    
    d(3,3,3)=1./8.*pow(1.+cb,3);
    d(3,3,2)=-sqrt(6.)/8.*sb*pow(1.+cb,2);
    d(3,3,1)=sqrt(15.)/8.*pow(sb,2)*(1.+cb);
    d(3,3,0)=-sqrt(5.)/4.*pow(sb,3);
    d(3,3,-1)=sqrt(15.)/8.*pow(sb,2)*(1-cb);
    d(3,3,-2)=-sqrt(6.)/8.*sb*pow(1-cb,2);
    d(3,3,-3)=1./8.*pow(1-cb,3);
    
    d(3,2,3)=pow(-1,2-3)*d(3,3,2);
    d(3,2,2)=-1./4.*pow(1+cb,2)*(2.-3.*cb);
    d(3,2,1)=sqrt(10.)/8.*sb*(1.-2.*cb-3.*pow(cb,2));
    d(3,2,0)=sqrt(30.)/4.*pow(sb,2)*cb;
    d(3,2,-1)=-sqrt(10.)/8.*sb*(1.+2.*cb-3*pow(cb,2));
    d(3,2,-2)=1./4.*pow(1-cb,2)*(2.+3.*cb);
    d(3,2,-3)=d(3,3,-2);

    d(3,1,3)=pow(-1,1-3)*d(3,3,1);
    d(3,1,2)=pow(-1,1-2)*d(3,2,1);
    d(3,1,1)=-1./8.*(1.+cb)*(1.+10.*cb-15.*pow(cb,2));
    d(3,1,0)=sqrt(3.)/4.*sb*(1.-5.*pow(cb,2));
    d(3,1,-1)=-1./8.*(1.-cb)*(1.-10.*cb-15.*pow(cb,2));
    d(3,1,-2)=d(3,2,-1);
    d(3,1,-3)=d(3,3,-1);

    d(3,0,3)=pow(-1,0-3)*d(3,3,0);
    d(3,0,2)=pow(-1,0-2)*d(3,2,0);
    d(3,0,1)=pow(-1,0-1)*d(3,1,0);
    d(3,0,0)=-1./2.*cb*(3.-5.*pow(cb,2));
    d(3,0,-1)=d(3,1,0);
    d(3,0,-2)=d(3,2,0);
    d(3,0,-3)=d(3,3,0);

    d(3,-1,3)=pow(-1,-1-3)*d(3,3,-1);
    d(3,-1,2)=pow(-1,-1-2)*d(3,2,-1);;
    d(3,-1,1)=pow(-1,-1-1)*d(3,1,-1);;
    d(3,-1,0)=d(3,0,1);
    d(3,-1,-1)=d(3,1,1);
    d(3,-1,-2)=d(3,2,1);
    d(3,-1,-3)=d(3,3,1);

    d(3,-2,3)=pow(-1,-2-3)*d(3,3,-2);;
    d(3,-2,2)=pow(-1,-2-2)*d(3,2,-2);;
    d(3,-2,1)=d(3,-1,2);
    d(3,-2,0)=d(3,0,2);
    d(3,-2,-1)=d(3,1,2);
    d(3,-2,-2)=d(3,2,2);
    d(3,-2,-3)=d(3,3,2);

    d(3,-3,3)=pow(-1,-3-3)*d(3,3,-3);
    d(3,-3,2)=pow(-1,-3-2)*d(3,2,-3);
    d(3,-3,1)=d(3,-1,3);
    d(3,-3,0)=d(3,0,3);
    d(3,-3,-1)=d(3,1,3);
    d(3,-3,-2)=d(3,2,3);
    d(3,-3,-3)=d(3,3,3);

    return d;
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_real::wigner_d(double beta)
{
    return wigner_d(beta,l_max);
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_real::wigner_d(double beta, int L)
{
    // M. Blanco and M. Florez and M Bermejo, "Evaluation of the rotation matrices in the basis of real spherical harmonics," Journal of Molecular Structure, 419, pp 19-27, 1997
    fdcl_FFTSO3_matrix<double> d(L);
    double cb, sb, sb2, cb2, tb2;
    int l, m, n;
    cb=cos(beta);
    sb=sin(beta);
    sb2=sin(beta/2.);
    cb2=cos(beta/2.);
    tb2=tan(beta/2.);
    
    d(0,0,0)=1.;
    
    d(1,0,0)=cb;
    d(1,1,-1)=pow(sb2,2.);  
    d(1,1,0)=-1./sqrt(2.)*sb;
    d(1,1,1)=pow(cb2,2.);
    
    // fill the lower triangular region
    for(l=2;l<=L;l++)
    {
        for(m=0;m<=l-2;m++)
        {
            for(n=-m;n<=m;n++)
            {
                d(l,m,n)=1./sqrt((pow(l,2.)-pow(m,2.))*(pow(l,2.)-pow(n,2.)))*
                    ( (l*(2.*l-1.)*d(1,0,0)-(2.*l-1.)*m*n/(l-1.))*d(l-1,m,n)
                    -sqrt((pow(l-1.,2.)-pow(m,2.))*(pow(l-1.,2)-pow(n,2.)))*l/(l-1.)*d(l-2,m,n) ); // (64)
            }
        }
        d(l,l,l)=d(1,1,1)*d(l-1,l-1,l-1); // (65)
        d(l,l-1,l-1)=(l*d(1,0,0)-l+1.)*d(l-1,l-1,l-1); //(66)
        
        for(n=l-1;n>=-l;n--)
            d(l,l,n)=-sqrt((l+n+1.)/(l-n))*tb2*d(l,l,n+1); // (67)
        
        for(n=l-2;n>=1-l;n--)
            d(l,l-1,n)=-(l*cb-n)/(l*cb-n-1.)*sqrt((l+n+1.)/(l-n))*tb2*d(l,l-1,n+1); // (68)
    }
    
    // fill remaining triangular regions with symmetry
    for(l=1;l<=L;l++)
    {   
        // upper triangle
        for(m=-l;m<0;m++)
            for(n=m;n<=-m;n++)
                d(l,m,n)=pow(-1.,m+n)*d(l,-m,-n);

        // left triangle
        for(n=-l;n<0;n++)
            for(m=n+1;m<-n;m++)
                d(l,m,n)=d(l,-n,-m);
        
        // right triangle
        for(n=1;n<=l;n++)
            for(m=-n+1;m<=n-1;m++)
                d(l,m,n)=pow(-1.,m+n)*d(l,n,m);
    }
    
    return d;
    
}

fdcl_FFTSO3_matrix_complex fdcl_FFTSO3_real::wigner_D(double alpha, double beta, double gamma, int L)
{
    fdcl_FFTSO3_matrix_real d(L);
    fdcl_FFTSO3_matrix_complex D(L);
    int l,m,n;
    
    d=wigner_d(beta,L);

    for(l=0;l<=L;l++)
        for(m=-l;m<=l;m++)
            for(n=-l;n<=l;n++)
                D(l,m,n)=d(l,m,n)*exp( -I*(alpha*((double)m) + gamma*((double)n)) );
    
    return D;
}

fdcl_FFTSO3_matrix_complex fdcl_FFTSO3_real::wigner_D_real_2(double alpha, double beta, double gamma, int L)
{
    fdcl_FFTSO3_matrix_complex D(L), C(L), D_real(L) ;
    int l;
    
    D=wigner_D(alpha,beta,gamma,L);
    C=matrix2rsph(L);

    for(l=0;l<=L;l++)
        D_real[l]=C[l].conjugate()*D[l]*C[l].transpose();

    return D_real;
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_real::wigner_D_real(double alpha, double beta, double gamma)
{
    return wigner_D_real(alpha,beta,gamma,l_max);
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_real::wigner_D_real(Matrix3 R)
{
    std::vector<double> abg;
    
    abg.resize(3);
    abg=R2Euler323(R);
    return wigner_D_real(abg[0],abg[1],abg[2]);
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_real::wigner_D_real_1(double alpha, double beta, double gamma, int L)
{
    fdcl_FFTSO3_matrix_real U(L), d(L);
    int l,m,n;
    std::vector<double> Phi;
    d=wigner_d(beta,L);

    Phi.resize(2);

    for(l=0;l<=L;l++)
    {
        for(m=-l;m<=l;m++)
        {
            for(n=-l;n<=l;n++)
            {
                Phi=compute_Phi(m,n,alpha,gamma);
                U(l,m,n)=Phi[0]*d(l,abs(m),abs(n))+Phi[1]*d(l,abs(m),-abs(n));
            }
        }
    }
       
    return U;
}

int fdcl_FFTSO3_real::signum(int x)
{
    int y;
    if (x > 0)
        y=1;
    else if (x < 0)
        y=-1;
    else
        y=0;

    return y;
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_real::wigner_D_real(double alpha, double beta, double gamma, int L)
{
    fdcl_FFTSO3_matrix_real U(L), d(L);
    int l,m,n;
    double cos_mamg, cos_ma_mg, sin_mamg, sin_ma_mg;    
    d=wigner_d(beta,L);

    for(l=0;l<=L;l++)
    {
        U(l,0,0)=d(l,0,0);
        for(m=1;m<=l;m++)
        {
            // U(l,-m,0)=pow(-1.,m)*sqrt(2.)*d(l,m,0)*sin(((double) m)*alpha);
            U(l,-m,0)=sqrt(2.)*d(l,-m,0)*sin(((double) m)*alpha);
            U(l,m,0)=pow(-1.,m)*sqrt(2.)*d(l,m,0)*cos(((double) m)*alpha);
            for(n=1;n<=l;n++)
            {
                cos_mamg = cos( ((double)m)*alpha + ((double)n)*gamma);
                cos_ma_mg = cos( ((double)m)*alpha - ((double)n)*gamma);
                sin_mamg = sin( ((double)m)*alpha + ((double)n)*gamma);
                sin_ma_mg = sin( ((double)m)*alpha - ((double)n)*gamma);

                U(l,m,n)=pow(-1.,m+n)*d(l,m,n)*cos_mamg + pow(-1.,m)*d(l,m,-n)*cos_ma_mg;
                U(l,m,-n)=pow(-1.,m)*d(l,m,-n)*sin_ma_mg
                    - pow(-1.,m+n)*d(l,m,n)*sin_mamg;
                U(l,-m,n)=-pow(-1.,n)*d(l,-m,n)*-sin_ma_mg
                    + d(l,-m,-n)*sin_mamg;
                U(l,-m,-n)=d(l,-m,-n)*cos_mamg
                    - pow(-1.,n)*d(l,-m,n)*cos_ma_mg;
            }
        }
        for(n=1;n<=l;n++)
        {
            // U(l,0,-n)=pow(-1.,n+1)*sqrt(2.)*d(l,0,n)*sin(((double) n)*gamma);
            U(l,0,-n)=-sqrt(2.)*d(l,0,-n)*sin(((double) n)*gamma);
            U(l,0,n)=pow(-1.,n)*sqrt(2.)*d(l,0,n)*cos(((double) n)*gamma);
        }
    }
        
    return U;
}

fdcl_FFTSO3_matrix_complex fdcl_FFTSO3_real::matrix2rsph(int L)
{
    fdcl_FFTSO3_matrix_complex C(L);
    int l,m;

    for(l=0;l<=L;l++)
    {
        // for(m=-l;m<0;m++)
        // {
            // C(l,m+l,m)=I/sqrt(2.);
            // C(l,m+l,-m)=-I/sqrt(2.)*pow(-1.,m);
        // }
        // C(l,0+l,0)=1.;
        // for(m=1;m<=l;m++)
        // {
            // C(l,-m,-m)=1./sqrt(2.);
            // C(l,-m,m)=pow(-1.,m)/sqrt(2.);
        // } 
        for(m=-l;m<0;m++)
        {
            C(l,m,m)=I/sqrt(2.);
            C(l,m,-m)=-I/sqrt(2.)*pow(-1.,m);
        }
        C(l,0,0)=1.;
        for(m=1;m<=l;m++)
        {
            C(l,m,-m)=1./sqrt(2.);
            C(l,m,m)=pow(-1.,m)/sqrt(2.);
        }
    }

    return C;
}

std::vector<double> fdcl_FFTSO3_real::compute_weight()
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
        
        sum*=1./((double)4*B*B*B)*sin((double)(2*j+1)*factor);
      
        weight[j]=sum;
    }
    
    return weight;
}

void fdcl_FFTSO3_real::check_weight()
{
    fdcl_FFTSO3_matrix<double> d(2*B-1);
    std::vector<double> sum;
    sum.resize(2*B);
    for (int l=0;l<2*B;l++)
        sum[l]=0.;

    this->compute_weight();
    
    for (int k=0;k<2*B; k++)
    {
        d=wigner_d(beta_k(k),2*B-1);
        for(int l=0;l<2*B;l++)
        {
            sum[l]+=d(l,0,0)*weight[k];
        }
    }
    
    cout << "fdcl_FFTSO3_real::check_weight" << endl;
    cout << "\\sum_k w_k d^l_00(beta_k) * 4B^2 = \\delta_{0,l}" << endl; 
    for (int l=0;l<2*B;l++)
        cout << "l=" << l << ": " << sum[l]*((double)4*B*B) << endl;
    
    
    int j1, j2, k, l;
    fdcl_FFTSO3_matrix_complex Delta(2*B-1);
    Delta.setZero();
    for(k=0;k<2*B;k++)
        for(j1=0;j1<2*B;j1++)
            for(j2=0;j2<2*B;j2++)
                for(l=0;l<2*B;l++)
                    Delta[l]+=weight[k]*wigner_D(alpha_j(j1),beta_k(k),gamma_j(j2),2*B-1)[l];
    
    
    cout << "\\sum_{j1,k,j2} w_k D(alpha_j1, beta_k, gamma_j2) = \\delta_{l,0}\\delta_{m,0}\\delta_{n,0}" << endl;
    for (int l=0;l<2*B;l++)
        cout << "l=" << l << ": " << Delta[l].norm() << endl;

}

void fdcl_FFTSO3_real::check_wigner_d()
{
    int l=3, N=1000;
    double beta;
    MatrixXd I, d_i_0, d_i_1;
    fdcl_FFTSO3_matrix_real d_beta(l), d_beta_explicit(l);
    
    I.resize(2*l+1,2*l+1);
    I.setIdentity();
    d_i_1.resize(2*l+1,2*l+1);
    d_i_1.setZero();
    d_i_0.resize(2*l-1,2*l-1);
    d_i_0.setZero();

    beta=(double)rand()/RAND_MAX;
    
    d_beta=wigner_d(beta,l);

    cout << "fdcl_FFTSO3_real::check_wigner_d" << endl;
    cout << "matrix orthogonality error: " << (d_beta[l].transpose()*d_beta[l]-I).norm() << endl;
    
    for(int i=0;i<N;i++)
    {
        beta=M_PI/((double)N)*((double)i);
        d_beta=wigner_d(beta,l);
        d_i_1+=d_beta[l].cwiseProduct(d_beta[l])*sin(beta)*M_PI/((double)N);
        d_i_0+=d_beta[l-1].cwiseProduct(d_beta[l].block(1,1,2*l-1,2*l-1))*sin(beta)*M_PI/((double)N);

    }
    d_i_1*=((double)(2*l+1))/2.;

    cout << "functional orthogonality error: " << endl;
    cout << d_i_0 << endl;
    cout << d_i_1 << endl;

    cout << "difference from the explicit expression" << endl;
    beta=(double)rand()/RAND_MAX*M_PI;
    
    d_beta=wigner_d(beta,l);
    d_beta_explicit=wigner_d_explicit(beta);
    
    cout << "beta =" << beta << endl;
    cout << "l=0: " << (d_beta[0]-d_beta_explicit[0]).norm() << endl;
    cout << "l=1: " << (d_beta[1]-d_beta_explicit[1]).norm() << endl;
    cout << "l=2: " << (d_beta[2]-d_beta_explicit[2]).norm() << endl;
    cout << "l=3: " << (d_beta[3]-d_beta_explicit[3]).norm() << endl;

    beta=-(double)rand()/RAND_MAX*M_PI;
    
    d_beta=wigner_d(beta,l);
    d_beta_explicit=wigner_d_explicit(beta);
    
    cout << "beta =" << beta << endl;
    cout << "l=0: " << (d_beta[0]-d_beta_explicit[0]).norm() << endl;
    cout << "l=1: " << (d_beta[1]-d_beta_explicit[1]).norm() << endl;
    cout << "l=2: " << (d_beta[2]-d_beta_explicit[2]).norm() << endl;
    cout << "l=3: " << (d_beta[3]-d_beta_explicit[3]).norm() << endl;

    cout << endl;
    
}

double fdcl_FFTSO3_real::delta(int i, int j)
{
    double delta=0.;
    
    if (i==j)
        delta=1.;
    
    return delta;
    
}

double fdcl_FFTSO3_real::inverse_transform(fdcl_FFTSO3_matrix_real F, double alpha, double beta, double gamma)
{
    double f=0.;
    int l,m,n;
    fdcl_FFTSO3_matrix_real U;
    
    U=wigner_D_real(alpha,beta,gamma);   
    
    for(l=0;l<=l_max;l++)
        for(m=-l;m<=l;m++)
            for(n=-l;n<=l;n++)
                f+=((double) 2*l+1 )* F(l,m,n)*U(l,m,n);
    
    return f;
}

double fdcl_FFTSO3_real::inverse_transform(fdcl_FFTSO3_matrix_real F, Matrix3 R)
{
    std::vector<double> abg;
    
    abg.resize(3);
    abg=R2Euler323(R);
    
    return inverse_transform(F,abg[0],abg[1],abg[2]); 
}

double fdcl_FFTSO3_real::beta_k(int k)
{
    return ((double)(2*k+1))*M_PI/4./((double)B);
}

double fdcl_FFTSO3_real::alpha_j(int j)
{
    return ((double)j)*M_PI/((double)B);
}

double fdcl_FFTSO3_real::gamma_j(int j)
{
    return alpha_j(j);
}

double fdcl_FFTSO3_real::f_4_check_forward_transform(double alpha, double beta, double gamma)
{
    double y=0.;
    int L=3;
    fdcl_FFTSO3_real tmp(L); 
    fdcl_FFTSO3_matrix_real U(L);
 
    U=tmp.wigner_D_real(alpha,beta,gamma);

    for(int l=0;l<=L;l++)
        for(int m=-l;m<=l;m++)
            for(int n=-l;n<=l;n++)
                y+=((double)2*l+1)*U(l,m,n);

    return y;
}

void fdcl_FFTSO3_real::check_forward_transform()
{
    fdcl_FFTSO3_matrix_real F(l_max), F_0(l_max);
    cout << "fdcl_FFTSO3_real::check_forward_transform_real" << endl;

    F=forward_transform(f_4_check_forward_transform);
    F_0=forward_transform_0(f_4_check_forward_transform);
    // 
    cout << "forward transform of a scaled real wigner D matrices: all of the Fourier parameters must be one for l \\leq 3, and zero otherwise." << endl;
// 
    cout << F ;
    cout << "error from forward_transform_0: " << (F-F_0).norm() << endl;
    cout << "fdcl_FFTSO3_real::check_forward_transform_real:completed" << endl << endl;
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_real::forward_transform(std::function <double(double, double, double)> func)
{
    fdcl_FFTSO3_matrix_real F_beta_1[2*B][2*B], F_beta_2[2*B][2*B], d_beta_k(l_max);
    fdcl_FFTSO3_matrix_real F(l_max);
    int j1, j2, k, l, m, n;
    std::vector<double> Phi;
    double f_j1kj2, alpha, beta, gamma;

    Phi.resize(2);
    compute_weight();

    for(j1=0;j1<2*B;j1++)
    {
        for(j2=0;j2<2*B;j2++)
        {
            F_beta_1[j1][j2].init(l_max);
            F_beta_2[j1][j2].init(l_max);
            F_beta_1[j1][j2].setZero();
            F_beta_2[j1][j2].setZero();
        }
    }

    for(k=0;k<2*B;k++)
    {
        beta=beta_k(k);     
        d_beta_k=wigner_d(beta_k(k));
        for(j1=0;j1<2*B;j1++)
        {
            alpha=alpha_j(j1);              
            for(j2=0;j2<2*B;j2++)
            {
                gamma=gamma_j(j2);
                f_j1kj2=func(alpha,beta,gamma);
                for (l=0;l<=l_max;l++)
                    for(m=-l;m<=l;m++)
                        for(n=-l;n<=l;n++)
                        {
                            F_beta_1[j1][j2](l,m,n)+=weight[k]*d_beta_k(l,abs(m),abs(n))*f_j1kj2; 
                            F_beta_2[j1][j2](l,m,n)+=weight[k]*d_beta_k(l,abs(m),-abs(n))*f_j1kj2; 
                        }
            }
        }
    }   

    F.setZero();
    for(j1=0;j1<2*B;j1++)
    {
        alpha=alpha_j(j1);              
        for(j2=0;j2<2*B;j2++)
        {
            gamma=gamma_j(j2);
            for (l=0;l<=l_max;l++)
                for(m=-l;m<=l;m++)
                    for(n=-l;n<=l;n++)
                    {
                        Phi=compute_Phi(m,n,alpha,gamma);
                        F(l,m,n)+=Phi[0]*F_beta_1[j1][j2](l,m,n)+Phi[1]*F_beta_2[j1][j2](l,m,n);
                    }
        }
    }
    return F;
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_real::forward_transform(std::function <double(Matrix3)> func) 
{
    // See lambda expression : https://www.geeksforgeeks.org/lambda-expression-in-c/
    return forward_transform([=] (double a, double b, double g)
            {
                return func(Euler3232R(a,b,g));
            });
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_real::forward_transform_0(std::function <double(double, double, double)> func)
{
    fdcl_FFTSO3_matrix_real F_beta_theta[2*B][2*B], F_beta_psi[2*B][2*B];
    fdcl_FFTSO3_matrix_real F_gamma_theta[2*B], F_gamma_psi[2*B];
    fdcl_FFTSO3_matrix_real F(l_max);
    int j1, j2, k, l, m, n;
    std::vector<fdcl_FFTSO3_matrix_real> TP;
    double f_j1kj2, alpha, beta, gamma;
    double cos_ng, sin_ng, sin_ma, cos_ma;

    TP.resize(2);
    TP[0].init(l_max);
    TP[1].init(l_max);

    compute_weight();

    for(j1=0;j1<2*B;j1++)
    {
        for(j2=0;j2<2*B;j2++)
        {
            F_beta_theta[j1][j2].init(l_max);
            F_beta_psi[j1][j2].init(l_max);
            F_beta_theta[j1][j2].setZero();
            F_beta_psi[j1][j2].setZero();
        }
    }

    for(k=0;k<2*B;k++)
    {
        beta=beta_k(k);     
        TP=compute_Theta_Psi(beta,l_max);
        for(j1=0;j1<2*B;j1++)
        {
            alpha=alpha_j(j1);              
            for(j2=0;j2<2*B;j2++)
            {
                gamma=gamma_j(j2);
                f_j1kj2=func(alpha,beta,gamma);
                for (l=0;l<=l_max;l++)
                    for(m=-l;m<=l;m++)
                        for(n=-l;n<=l;n++)
                        {
                            F_beta_theta[j1][j2](l,m,n)+=weight[k]*TP[0](l,m,n)*f_j1kj2; 
                            F_beta_psi[j1][j2](l,m,n)+=weight[k]*TP[1](l,m,n)*f_j1kj2; 
                        }
            }
        }
    }   

    for(j1=0;j1<2*B;j1++)
    {
        F_gamma_theta[j1].init(l_max);
        F_gamma_psi[j1].init(l_max);
        F_gamma_theta[j1].setZero();
        F_gamma_psi[j1].setZero();
    }

    for(j1=0;j1<2*B;j1++)
        for(j2=0;j2<2*B;j2++)
        {
            gamma=gamma_j(j2);
            for(l=0;l<=l_max;l++)
                for(m=-l;m<=l;m++)
                {
                    n=0;
                    if(m >= 0)
                    {
                        F_gamma_psi[j1](l,m,n)+=F_beta_psi[j1][j2](l,m,n);
                    }
                    else
                    {
                        F_gamma_theta[j1](l,m,n)+=F_beta_theta[j1][j2](l,m,n);
                    }


                    for(n=1;n<=l;n++)
                    {
                        cos_ng=cos( ((double)n)*gamma);
                        sin_ng=sin( ((double)n)*gamma);

                        if(m>=0)
                        {
                            F_gamma_theta[j1](l,m,n)+=sin_ng*F_beta_theta[j1][j2](l,m,n);
                            F_gamma_psi[j1](l,m,n)+=cos_ng*F_beta_psi[j1][j2](l,m,n);
                            F_gamma_theta[j1](l,m,-n)+=cos_ng*F_beta_theta[j1][j2](l,m,-n);
                            F_gamma_psi[j1](l,m,-n)+=-sin_ng*F_beta_psi[j1][j2](l,m,-n);
                        }
                        else
                        {
                            F_gamma_theta[j1](l,m,n)+=cos_ng*F_beta_theta[j1][j2](l,m,n);
                            F_gamma_psi[j1](l,m,n)+=sin_ng*F_beta_psi[j1][j2](l,m,n);
                            F_gamma_theta[j1](l,m,-n)+=-sin_ng*F_beta_theta[j1][j2](l,m,-n);
                            F_gamma_psi[j1](l,m,-n)+=cos_ng*F_beta_psi[j1][j2](l,m,-n);
                        }

                        // if( (m>=0 && n >=0) || (m<0 && n<0))
                        // {
                            // F_gamma_theta[j1](l,m,n)+=sin_ng*F_beta_theta[j1][j2](l,m,n);
                            // F_gamma_psi[j1](l,m,n)+=cos_ng*F_beta_psi[j1][j2](l,m,n);
                        // }
                        // else
                        // {
                            // F_gamma_theta[j1](l,m,n)+=cos_ng*F_beta_theta[j1][j2](l,m,n);
                            // F_gamma_psi[j1](l,m,n)+=sin_ng*F_beta_psi[j1][j2](l,m,n);
                        // }
                    }
                }
        }

    F.setZero();
    for(j1=0;j1<2*B;j1++)
    {
        alpha=alpha_j(j1);              
        for (l=0;l<=l_max;l++)
        {
            m=0;
            for(n=-l;n<=l;n++)
                F(l,m,n)+=F_gamma_psi[j1](l,m,n);

            for(m=1;m<=l;m++)
            {
                sin_ma=sin( ((double) m)*alpha );
                cos_ma=cos( ((double) m)*alpha );

                for(n=-l;n<=l;n++)
                {
                    F(l,m,n)+=sin_ma*F_gamma_theta[j1](l,m,n)+cos_ma*F_gamma_psi[j1](l,m,n);
                    F(l,-m,n)+=-sin_ma*F_gamma_theta[j1](l,-m,n)+cos_ma*F_gamma_psi[j1](l,-m,n);
                }
            }
        }
            // }
            // 
            // for(m=-l;m<=l;m++)
            // {
                // sin_ma=sin( ((double) m)*alpha );
                // cos_ma=cos( ((double) m)*alpha );
// 
                // for(n=-l;n<=l;n++)
                    // F(l,m,n)+=sin_ma*F_gamma_theta[j1](l,m,n)+cos_ma*F_gamma_psi[j1](l,m,n);
            // }
    }

    return F;
}

std::vector<double> fdcl_FFTSO3_real::compute_Phi(int m, int n, double alpha, double gamma)
{
    std::vector<double> Phi;
    Phi.resize(2);

    if (m*n > 0)
    {
        Phi[0] = pow(-1.,m-n) * cos( ((double)m)*alpha + ((double)n)*gamma);
        Phi[1] = pow(-1.,m)* signum(m) *  cos( ((double)m)*alpha - ((double)n)*gamma);
    }
    else if (m*n < 0)
    {
        Phi[0] = -pow(-1.,m-n)*sin( ((double)m)*alpha - ((double)n)*gamma);
        Phi[1] = pow(-1.,m)* signum(m) *sin( ((double)m)*alpha + ((double)n)*gamma);
    }
    else if ( (m>0 && n==0) || (m==0 && n > 0) )
    {
        Phi[0] = pow(-1.,m-n)*sqrt(2)* cos( ((double)m)*alpha + ((double)n)*gamma);
        Phi[1] = 0.0;
    }
    else if ( (m<0 && n==0) || (m==0 & n < 0) )
    {
        Phi[0] = -pow(-1.,m-n)*sqrt(2.)*sin( ((double)m)*alpha - ((double)n)*gamma);
        Phi[1] = 0.0;
    }
    else if(m==0 && n==0)
    {
        Phi[0] = 1.0;
        Phi[1] = 0.0;
    }  

    return Phi;
}

std::vector<fdcl_FFTSO3_matrix_real> fdcl_FFTSO3_real::compute_Theta_Psi(double beta, int L)
{
    std::vector<fdcl_FFTSO3_matrix_real> TP;
    fdcl_FFTSO3_matrix_real d(L);
    int l,m,n;
    double A, B, C;

    TP.resize(2);
    TP[0].init(L);
    TP[1].init(L);

    d=wigner_d(beta,L);

    for(l=0;l<=L;l++)
    {
        TP[0](l,0,0)=0.;
        TP[1](l,0,0)=d(l,0,0);

        m=0;
        for(n=1;n<=l;n++)
        {
            C=sqrt(2.)*d(l,abs(m),abs(n));
            TP[0](l,m,n)=-pow(-1.,m-n)*C;
            TP[1](l,m,n)=pow(-1,m-n)*C;
            TP[0](l,m,-n)=TP[0](l,m,n);
            TP[1](l,m,-n)=TP[1](l,m,n);
        }
        for(m=1;m<=l;m++)
        {
            n=0;
            C=sqrt(2.)*d(l,abs(m),abs(n));
            TP[0](l,m,n)=-pow(-1.,m-n)*C;
            TP[1](l,m,n)=pow(-1.,m-n)*C;
            TP[0](l,-m,n)=TP[0](l,m,n);
            TP[1](l,-m,n)=TP[1](l,m,n);
            for(n=1;n<=l;n++)
            {
                A=pow(-1.,m-n)*d(l,abs(m),abs(n));
                B=pow(-1.,m)*signum(m)*d(l,abs(m),-abs(n));
                TP[0](l,m,n)=-A+B;
                TP[1](l,m,n)=A+B;
                TP[0](l,-m,n)=-A-B;
                TP[1](l,-m,n)=A-B;

                TP[0](l,m,-n)=TP[0](l,m,n);
                TP[1](l,m,-n)=TP[1](l,m,n);
                TP[0](l,-m,-n)=TP[0](l,-m,n);
                TP[1](l,-m,-n)=TP[1](l,-m,n);
            }
        }
    }


    return TP;
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_real::wigner_D_real_0(double alpha, double beta, double gamma, int L)
{
    fdcl_FFTSO3_matrix_real U(L);
    std::vector<fdcl_FFTSO3_matrix_real> TP(L);
    int l,m,n;
    double cos_ma, sin_ma, cos_ng, sin_ng;
    TP.resize(2);
    // TP[1]=Theta, TP[2]=Psi;
    TP=compute_Theta_Psi(beta,L);
    
    for(l=0;l<=L;l++)
    {
        U(l,0,0)=TP[1](l,0,0);
        for(n=1;n<=l;n++)
        {
            // when m=0
            sin_ng=sin( ((double)n)*gamma);
            cos_ng=cos( ((double)n)*gamma);
            U(l,0,n)=cos_ng*TP[1](l,0,n);
            U(l,0,-n)=-sin_ng*TP[1](l,0,-n);
        }
        for(m=1;m<=l;m++)
        {
            // when n=0
            cos_ma=cos( ((double)m)*alpha);
            sin_ma=sin( ((double)m)*alpha);
            U(l,m,0)=cos_ma*TP[1](l,m,0);
            U(l,-m,0)=-sin_ma*TP[0](l,-m,0);
            for(n=1;n<=l;n++)
            {
                sin_ng=sin( ((double)n)*gamma);
                cos_ng=cos( ((double)n)*gamma);
                U(l,m,n)=sin_ma*sin_ng*TP[0](l,m,n)+cos_ma*cos_ng*TP[1](l,m,n);
                U(l,-m,-n)=sin_ma*sin_ng*TP[0](l,-m,-n)+cos_ma*cos_ng*TP[1](l,-m,-n);
                U(l,m,-n)=sin_ma*cos_ng*TP[0](l,m,-n)-cos_ma*sin_ng*TP[1](l,m,-n);
                U(l,-m,n)=-sin_ma*cos_ng*TP[0](l,-m,n)+cos_ma*sin_ng*TP[1](l,-m,n);
            }   
        }
    }

    return U;
}

void fdcl_FFTSO3_real::check_wigner_D_real()
{
    int L=10;
    double alpha, beta, gamma;
    fdcl_FFTSO3_matrix_complex U_2(L);
    fdcl_FFTSO3_matrix_real U(L), U_0(L), U_1(L);
    fdcl_tictoc tictoc;

    alpha=(double)rand()/RAND_MAX*M_PI*2.;
    beta=(double)rand()/RAND_MAX*M_PI;
    gamma=(double)rand()/RAND_MAX*M_PI*2.;

    cout << "fdcl_FFTSO3_real::check_wigner_D_real" << endl;
    cout << "alpha=" << alpha << ", beta=" << beta << ", gamma=" << gamma << endl;
    U=wigner_D_real(alpha,beta,gamma,L);
    U_0=wigner_D_real_0(alpha,beta,gamma,L);
    U_1=wigner_D_real_1(alpha,beta,gamma,L);
    U_2=wigner_D_real_2(alpha,beta,gamma,L);

    cout << "error from wigner_D_real_0: " << (U-U_0).norm() << endl;
    cout << "error from wigner_D_real_1: " << (U-U_1).norm() << endl;
    cout << "error from wigner_D_real_2: " << (U-U_2.real()).norm() << endl;

    tictoc.tic();
    for(int i=0;i<=500;i++)
        wigner_D_real(alpha,beta,gamma,L);
    tictoc.toc("wigner_D_real");
    
    tictoc.tic();
    for(int i=0;i<=500;i++)
        wigner_D_real_0(alpha,beta,gamma,L);
    tictoc.toc("wigner_D_real_0");
    
    tictoc.tic();
    for(int i=0;i<=500;i++)
        wigner_D_real_1(alpha,beta,gamma,L);
    tictoc.toc("wigner_D_real_1");
    
    tictoc.tic();
    for(int i=0;i<=500;i++)
        wigner_D_real_2(alpha,beta,gamma,L);
    tictoc.toc("wigner_D_real_2");
    
    cout << "fdcl_FFTSO3_real::check_wigner_D_real completed" << endl << endl;
}

void fdcl_FFTSO3_real::check_Clebsch_Gordon()
{
    int l1=1, l2=4, l, m1, m2, n1, n2;
    double alpha, beta, gamma;
    fdcl_FFTSO3_matrix_real U(l1+l2);
    double y, y_CB=0.;
    double error=0.;
    std::vector<int> M, N;

    alpha=(double)rand()/RAND_MAX*2.*M_PI;
    beta=(double)rand()/RAND_MAX*M_PI;
    gamma=(double)rand()/RAND_MAX*2.*M_PI;

    U=wigner_D_real(alpha,beta,gamma,l1+l2);

    c.compute(l1,l2);
    double cr,ci;
    for(int i=0;i<(2*l1+1)*(2*l2+1);i++)
        for(int j=0;j<(2*l1+1)*(2*l2+1);j++)
        {
            if (abs(real(c.c(i,j))) < 1e-10)
                cr=0.;
            else
                cr=real(c.c(i,j));

            if (abs(imag(c.c(i,j))) < 1e-10)
                ci=0.;
            else
                ci=imag(c.c(i,j));

            c.c(i,j)=cr+I*ci;
        }




    cout << c.c << endl;

    cout << "fdcl_FFTSO3_real:check_Clebsch_Gordon" << endl;
    cout << "l1 = " << l1 << ", l2 = " << l2 << endl;
    cout << "alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << endl;

    for(m1=-l1;m1<=l1;m1++)
        for(n1=-l1;n1<=l1;n1++)
            for(m2=-l2;m2<=l2;m2++)
                for(n2=-l2;n2<=l2;n2++)
                {
                    M.clear();
                    M.insert(M.end(),{m1+m2,m1-m2,-m1+m2,-m1-m2});
                    std::sort(M.begin(),M.end());
                    auto last_M=std::unique(M.begin(),M.end());
                    M.erase(last_M,M.end());

                    N.clear();
                    N.insert(N.end(),{n1+n2,n1-n2,-n1+n2,-n1-n2});
                    std::sort(N.begin(),N.end());
                    auto last_N=std::unique(N.begin(),N.end());
                    N.erase(last_N,N.end());

                    y=U(l1,m1,n1)*U(l2,m2,n2);
                    y_CB=0.;
                    for(int m : M)
                        for(int n : N)
                            for(l=max(max(abs(l1-l2),abs(m)),abs(n));l<=l1+l2;l++)
                            {
                                y_CB+=real(c(l,m,l1,m1,l2,m2)*std::conj(c(l,n,l1,n1,l2,n2)))*U(l,m,n);

                                double tmp=std::imag(c(l,m,l1,m1,l2,m2)*std::conj(c(l,n,l1,n1,l2,n2))); 
                                if (abs(tmp) > 1e-6)
                                    cout << tmp << endl;
                            }

                    // cout << "y = " << y << endl;
                    // cout << "y_Clebsch_Gordon = " << y_CB << endl;
                    if(abs(y-y_CB) > error)
                        error = abs(y-y_CB);
                }

    cout << "error = " << error << endl;
    cout << "fdcl_FFTSO3_real:check_Clebsch_Gordon:completed" << endl << endl;
}


