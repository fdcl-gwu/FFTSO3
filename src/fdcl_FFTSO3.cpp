#include "fdcl_FFTSO3.hpp"

fdcl_FFTSO3_complex::fdcl_FFTSO3_complex(int l_max)
{
    init(l_max);
}

void fdcl_FFTSO3_complex::init(int l_max)
{
    this->l_max=l_max;  
    this->B=l_max+1;
    d_beta.resize(2*B);
    weight.resize(2*B);
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_complex::wigner_d_explicit(double beta)
{
    // D. Varshalovich, A. Moskalev, and V. Khersonskii, Quantum Theory of Angular Momentum, World Scientific, 1988, Chapter 4
    
    fdcl_FFTSO3_matrix<double> d(3);
    double cb, sb;
    cb=cos(beta);
    sb=sin(beta);
   
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

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_complex::wigner_d(double beta)
{
    return wigner_d(beta,l_max);
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_complex::wigner_d(double beta, int L)
{
    // M. Blanco and M. Florez and M Bermejo, "Evaluation of the rotation matrices in the basis of real spherical harmonics," Journal of Molecular Structure, 419, pp 19-27, 1997
    fdcl_FFTSO3_matrix<double> d(L);
    double cb, sb;
    double sb2, cb2, tb2;
    int l,m,n;

    cb=cos(beta);
    sb=sin(beta);
    sb2=sin(beta/2.);
    cb2=cos(beta/2.);
    tb2=tan(beta/2.);
   
    d(0,0,0)=1.;
   
	if(L>=1)
	{ 
		d(1,0,0)=cb;
		d(1,1,-1)=pow(sb2,2.);  
		d(1,1,0)=-1./sqrt(2.)*sb;
		d(1,1,1)=pow(cb2,2.);
	}
    
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

fdcl_FFTSO3_matrix_complex fdcl_FFTSO3_complex::wigner_D(double alpha, double beta, double gamma)
{
    return wigner_D(alpha,beta,gamma,l_max);
}

fdcl_FFTSO3_matrix_complex fdcl_FFTSO3_complex::wigner_D(double alpha, double beta, double gamma, int L)
{
    fdcl_FFTSO3_matrix_real d(L);
    fdcl_FFTSO3_matrix_complex D(L);
    std::complex<double> expIMA;
    int l,m,n;
    
    d=wigner_d(beta,L);

    for(l=0;l<=L;l++)
    {
#pragma omp parallel for private(expIMA,m,n)
        for(m=-l;m<=l;m++)
        {
            expIMA=exp(-I*(double)m*alpha);
            for(n=-l;n<=l;n++)
            {
                D(l,m,n)=d(l,m,n)*expIMA*exp( -I*(gamma*((double)n)) );
            }
        }
    }
    
    return D;
}

fdcl_FFTSO3_matrix_complex fdcl_FFTSO3_complex::wigner_D(Matrix3 R)
{
    std::vector<double> abg;
    
    abg.resize(3);
    abg=R2Euler323(R);
    
    return wigner_D(abg[0],abg[1],abg[2]);      
}

std::vector<double> fdcl_FFTSO3_complex::compute_weight()
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

double fdcl_FFTSO3_complex::check_weight()
{
    double error=1.;
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
    
    if(check_verbose)
    {
        cout << "fdcl_FFTSO3_complex::check_weight" << endl;
        cout << "\\sum_k w_k d^l_00(beta_k) * 4B^2 = \\delta_{0,l}" << endl; 
        for (int l=0;l<2*B;l++)
            cout << "l=" << l << ": " << sum[l]*((double)4*B*B) << endl;
    }

    error = abs(sum[0]*((double)4*B*B)-1.);
    for (int l=1; l<2*B; l++)
        error += abs(sum[l]); 
    
    int j1, j2, k, l;
    fdcl_FFTSO3_matrix_complex Delta(2*B-1);
    Delta.setZero();
    for(k=0;k<2*B;k++)
        for(j1=0;j1<2*B;j1++)
            for(j2=0;j2<2*B;j2++)
                for(l=0;l<2*B;l++)
                    Delta[l]+=weight[k]*wigner_D(alpha_j(j1),beta_k(k),gamma_j(j2),2*B-1)[l];
    
    if(check_verbose)
    {
        cout << "\\sum_{j1,k,j2} w_k D(alpha_j1, beta_k, gamma_j2) = \\delta_{l,0}\\delta_{m,0}\\delta_{n,0}" << endl;
        for (int l=0;l<2*B;l++)
            cout << "l=" << l << ": " << Delta[l].norm() << endl;
    }

    error += abs(Delta(0,0,0)-1.0);
    for (int l=1; l<2*B; l++)
        error += Delta[l].norm();
    
    cout << "fdcl_FFTSO3_complex::check_weight: error = " << error << endl;
    return error;

}

double fdcl_FFTSO3_complex::check_wigner_d()
{
    int L=3, N=2000;
    double beta;
    MatrixXd I, d_i_0, d_i_1;
    fdcl_FFTSO3_matrix_real d_beta(L), d_beta_explicit(L);
    Eigen::VectorXd error(4);
    
    error.setZero();
    I.resize(2*L+1,2*L+1);
    I.setIdentity();
    d_i_1.resize(2*L+1,2*L+1);
    d_i_1.setZero();
    d_i_0.resize(2*L-1,2*L-1);
    d_i_0.setZero();

    beta=(double)rand()/RAND_MAX;
    
    d_beta=wigner_d(beta,L);

    error(0) = (d_beta[L].transpose()*d_beta[L]-I).norm();
    if(check_verbose)
    {
        cout << "fdcl_FFTSO3_complex::check_wigner_d" << endl;
        cout << "matrix orthogonality error: " << error << endl;
    }
   

    if(check_verbose)
        cout << "functional orthogonality error: " << endl;

    for(int i=0;i<N;i++)
    {
        beta=M_PI/((double)N)*((double)i);
        d_beta=wigner_d(beta,L);
        d_i_1+=d_beta[L].cwiseProduct(d_beta[L])*sin(beta)*M_PI/((double)N);
        d_i_0+=d_beta[L-1].cwiseProduct(d_beta[L].block(1,1,2*L-1,2*L-1))*sin(beta)*M_PI/((double)N);
    }
    d_i_1*=((double)(2*L+1))/2.;

    if(check_verbose)
    {
        cout << d_i_0 << endl;
        cout << d_i_1 << endl;
    }
    error(1) = d_i_0.norm();
    d_i_1.array() -=1.0;
    error(1) += d_i_1.norm();


    if(check_verbose)
        cout << "difference from the explicit expression" << endl;

    beta=(double)rand()/RAND_MAX*M_PI;
    
    d_beta=wigner_d(beta,L);
    d_beta_explicit=wigner_d_explicit(beta);
     
    if(check_verbose)
    {
        cout << "beta =" << beta << endl;
        for (int l=0; l<=3; l++)
            cout << "l=" << l << ": " << (d_beta[l]-d_beta_explicit[l]).norm() << endl;
    }
    
    for (int l=0; l<=3; l++)
        error(2) += (d_beta[l]-d_beta_explicit[l]).norm();
    

    beta=-(double)rand()/RAND_MAX*M_PI;
    
    d_beta=wigner_d(beta,L);
    d_beta_explicit=wigner_d_explicit(beta);
    
    if(check_verbose)
    {
        cout << "beta =" << beta << endl;
        for (int l=0; l<=3; l++)
            cout << "l=" << l << ": " << (d_beta[l]-d_beta_explicit[l]).norm() << endl;
    }
    for (int l=0; l<=3; l++)
        error(3) += (d_beta[l]-d_beta_explicit[l]).norm();

    cout << "fdcl_FFTSO3_complex::check_wigner_d: error = " << error.maxCoeff()*1.e-9 << endl;
    
    return error.maxCoeff()*1.e-9;
   
}

std::vector<double> fdcl_FFTSO3_complex::character(double theta)
{
    std::vector<double> chi;
    fdcl_FFTSO3_matrix_real d(l_max);
        
    chi.resize(l_max+1);
    
    for(int l=0;l<=l_max;l++)
        chi[l]=sin( ((double)2*l+1)/2.*theta )/sin(theta/2.);
    
    return chi;
} 
    
std::vector<fdcl_FFTSO3_matrix_complex> fdcl_FFTSO3_complex::deriv_wigner_D()
{
    std::vector<fdcl_FFTSO3_matrix_complex> u;
    double c_n, cn;
    u.resize(4);
    u[1].init(l_max);
    u[2].init(l_max);
    u[3].init(l_max);   
    
    int l,m,n;
    
    for(l=0;l<=l_max;l++)
    {
        for(m=-l;m<=l;m++)
        {
            n=m;
            u[3](l,m,n)=-I*((double)m);
        }
        for(m=-l;m<l;m++)
        {
            n=m+1;
            c_n=sqrt(((double)l+n)*((double)l-n+1));            
            u[2](l,m,n)=0.5*c_n;
            u[1](l,m,n)=-0.5*I*c_n;
            
        }
        for(m=-l+1;m<=l;m++)
        {
            n=m-1;
            cn=sqrt(((double)l-n)*((double)l+n+1));         
            u[2](l,m,n)=-0.5*cn;
            u[1](l,m,n)=-0.5*I*cn;                      
        }
    }
    
    return u;
}

double fdcl_FFTSO3_complex::check_deriv_wigner_D()
{
    std::vector<fdcl_FFTSO3_matrix_complex> u;
    fdcl_FFTSO3_matrix_complex D, D_new;
    std::vector<double>abg;
    Eigen::Matrix<double, 3, 1> ei;
    double eps=1.e-6;
    Eigen::VectorXd error(l_max);
    error.setZero();
    
    u=deriv_wigner_D();
    D=wigner_D(0,0,0);

    if(check_verbose)
        cout << "fdcl_FFTSO3_complex::check_deriv_wigner_D";

    for(int i=1;i<=3;i++)
    {
        if(check_verbose)
            cout << endl << "u_" << i << endl;

        ei.setZero();
        ei(i-1)=1.;
        abg=R2Euler323(expm_SO3(ei*eps));
        D_new=wigner_D(abg[0],abg[1],abg[2]);

        for(int l=1;l<=l_max;l++)
        {
            if(check_verbose)
                cout << "l=" << l << ", error= " << ((D_new[l]-D[l])/eps-u[i][l]).norm()/u[i][l].norm() <<  endl;

            error(l-1)+=((D_new[l]-D[l])/eps-u[i][l]).norm()/u[i][l].norm();
//          cout << u[i][l] << endl << endl; // analytic derivative
//          cout << (D_new[l]-D[l])/eps << endl; // numerical derivative
        }
    }

    cout << "fdcl_FFTSO3_complex::check_deriv_wigner_D: error = " << error.maxCoeff()*1.e-9 << endl;

    return error.maxCoeff()*1.e-9;

}

complex<double> fdcl_FFTSO3_complex::inverse_transform(fdcl_FFTSO3_matrix_complex F, double alpha, double beta, double gamma)
{
	init(F.l_max);
	fdcl_FFTSO3_matrix_real d(l_max);
	Eigen::MatrixXcd Fmn(2*B,2*B);
	Eigen::VectorXcd expIMA(2*B), expING(2*B);

    d=wigner_d(beta);
    Fmn.setZero();

#pragma omp parallel
    {
#pragma omp for
        for(int m=-l_max; m<=l_max; m++)
        {
            expIMA(m+l_max)=exp(-I*(double)m*alpha);
            expING(m+l_max)=exp(-I*(double)m*gamma);
        }

        Eigen::MatrixXcd Fmn_local(2*B,2*B);
        Fmn_local.setZero();
#pragma omp for
        for(int m=-l_max; m<=l_max; m++)
            for(int n=-l_max; n<=l_max; n++)
            {
                for(int l=max(abs(m),abs(n)); l<=l_max; l++)
                    Fmn_local(m+l_max,n+l_max)+=(double)(2*l+1)*F(l,m,n)*d(l,m,n);
            }
#pragma omp critical
        Fmn=Fmn+Fmn_local;
    }

    return expIMA.transpose()*Fmn*expING;

	// alternative method: direct sum: slow
    // complex<double> f=0.;
    // int l,m,n;
    // init(F.l_max);
    // fdcl_FFTSO3_matrix_complex D;
    // 
    // D=wigner_D(alpha,beta,gamma);   
    // 
    // for(l=0;l<=l_max;l++)
        // for(m=-l;m<=l;m++)
            // for(n=-l;n<=l;n++)
                // f+=((double) 2*l+1 )* F(l,m,n)*D(l,m,n);
    // 
    // return f;

	// alternative method: smaller memory requirement, but little slower
	// complex<double> f=0.;

	// int l,m,n;
	// init(F.l_max);
	// fdcl_FFTSO3_matrix_real d(l_max);
	// complex<double> tmp={0.,0.};
// 
	// d=wigner_d(beta);
    // for(m=-l_max; m<=l_max; m++)
		// for(n=-l_max; n<=l_max; n++)
		// {
			// tmp=0.;
   			// for(l=max(abs(m),abs(n)); l<=l_max; l++)
				// tmp+=(double)(2*l+1)*F(l,m,n)*d(l,m,n);
			// f+=exp(-I*((double)m*alpha+(double)n*gamma))*tmp;
		// }
// 
    // return f;
}

complex<double> fdcl_FFTSO3_complex::inverse_transform(fdcl_FFTSO3_matrix_complex F, Matrix3 R)
{
    std::vector<double> abg;
    
    abg.resize(3);
    abg=R2Euler323(R);
    
    return inverse_transform(F,abg[0],abg[1],abg[2]);
}

double fdcl_FFTSO3_complex::beta_k(int k)
{
    return ((double)(2*k+1))*M_PI/4./((double)B);
}

double fdcl_FFTSO3_complex::alpha_j(int j)
{
    return ((double)j)*M_PI/((double)B);
}

double fdcl_FFTSO3_complex::gamma_j(int j)
{
    return alpha_j(j);
}

complex<double> fdcl_FFTSO3_complex::f_4_check_transform(double alpha, double beta, double gamma)
{
    return inverse_transform(F_4_check,alpha,beta,gamma);
}

double fdcl_FFTSO3_complex::check_transform()
{
    double error=1.;
    F_4_check.init(l_max);
    F_4_check.setRandom();

    auto func= std::bind(&fdcl_FFTSO3_complex::f_4_check_transform, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);

    error = (F_4_check-forward_transform(func)).norm();
    cout << "fdcl_FFTSO3_complex::check_transform: l_max=" << l_max << " : error = " << error << endl;

    return error;
}

fdcl_FFTSO3_matrix_complex fdcl_FFTSO3_complex::forward_transform_0(std::function <complex<double>(double, double, double)> func)
{
    int j1, j2, k, l, m, n;
    fdcl_FFTSO3_matrix_real d_beta_k(l_max);
    fdcl_FFTSO3_matrix_complex F(l_max);
    
    complex<double> exp_imalpha, exp_ingamma, f_j1kj2;
    double alpha, beta, gamma;
    compute_weight();
    
    F.setZero();
   
    for(k=0;k<2*B;k++)  
    {
        beta=beta_k(k);
        d_beta_k=wigner_d(beta);
        for(j1=0;j1<2*B;j1++)
        {
            alpha=alpha_j(j1);
            for(j2=0;j2<2*B;j2++)
            {
                gamma=gamma_j(j2);
                f_j1kj2=func(alpha,beta,gamma);
                for(l=0;l<B;l++)
                {
                    for(m=-l;m<=l;m++)
                    {
                        exp_imalpha=exp(I*((double)m)*alpha);
                        for(n=-l;n<=l;n++)
                        {
                            exp_ingamma=exp(I*((double)n)*gamma);
                            F(l,m,n)+=weight[k]*exp_imalpha*f_j1kj2*exp_ingamma*d_beta_k(l,m,n);
                        }
                    }
                }
            }
        }
    }
        
    return F;
}

fdcl_FFTSO3_matrix_complex fdcl_FFTSO3_complex::forward_transform(std::function <complex<double>(double, double, double)> func)
{
    return this->forward_transform(func,0);
}

fdcl_FFTSO3_matrix_complex fdcl_FFTSO3_complex::forward_transform(std::function <complex<double>(double, double, double)> func, bool is_real)
{
    fdcl_FFTSO3_matrix_complex F(l_max);
	F.setZero();
	compute_weight();

#pragma omp parallel 
    {
        fdcl_FFTSO3_matrix_complex F_local(l_max);
        fdcl_FFTSO3_matrix_real d_beta_k(l_max);
        Eigen::MatrixXcd func_k(2*B,2*B), F_k(2*B,2*B), F_k_new(2*B,2*B);
        Eigen::VectorXcd tmp_out(2*B);
        int j1, j2, k, l, m, n;
        double alpha, beta;
        Eigen::FFT<double> fft;

#pragma omp for 
        for(k=0;k<2*B;k++)
        {
            beta=beta_k(k);     

            for(j1=0;j1<2*B;j1++)
            {
                alpha=alpha_j(j1);              
                for(j2=0;j2<2*B;j2++)
                    func_k(j1,j2)=func(alpha,beta,gamma_j(j2));
            }	

            F_k.setZero();

            for(j1=0;j1<2*B;j1++)
            {
                fft.fwd(tmp_out, func_k.row(j1));
                F_k.row(j1)=tmp_out;
            }
            for(j2=0;j2<2*B;j2++)
            {
                fft.fwd(tmp_out, F_k.col(j2));
                F_k.col(j2)=tmp_out.transpose();
            }

            d_beta_k=wigner_d(beta,l_max);
            
            F_local.setZero();

            if(!is_real)
            {
                for(l=0; l<=l_max; l++)
                {
                    for(m=1; m<=l; m++)
                    {
                        for(n=1; n<=l; n++)	
                            F_local(l,m,n)+=weight[k]*d_beta_k(l,m,n)*F_k(2*B-m,2*B-n);
                        for(n=-l; n<=0; n++)	
                            F_local(l,m,n)+=weight[k]*d_beta_k(l,m,n)*F_k(2*B-m,-n);
                    }
                    for(m=-l; m<=0; m++)
                    {
                        for(n=1; n<=l; n++)	
                            F_local(l,m,n)+=weight[k]*d_beta_k(l,m,n)*F_k(-m,2*B-n);
                        for(n=-l; n<=0; n++)	
                            F_local(l,m,n)+=weight[k]*d_beta_k(l,m,n)*F_k(-m,-n);
                    }
                }
            }
#pragma omp critical
            F=F+F_local;
        }
    }

	return F;
}

fdcl_FFTSO3_matrix_complex fdcl_FFTSO3_complex::forward_transform_1(std::function <complex<double>(double, double, double)> func)
{
    fdcl_FFTSO3_matrix_complex F_beta[2*B][2*B];
    int j1, j2, k, l, m, n;
    fdcl_FFTSO3_matrix_real d_beta_k(l_max);
    fdcl_FFTSO3_matrix_complex F_gamma[2*B];
    fdcl_FFTSO3_matrix_complex F(l_max);
    
    complex<double> exp_imalpha, exp_ingamma, f_j1kj2;
    double alpha, beta, gamma;
    
    compute_weight();
    
    for(j1=0;j1<2*B;j1++)
    {
        for(j2=0;j2<2*B;j2++)
        {
            F_beta[j1][j2].init(l_max);
            F_beta[j1][j2].setZero();
        }
    }

    for(k=0;k<2*B;k++)
    {
        beta=beta_k(k);     
        d_beta_k=wigner_d(beta_k(k),l_max);
        for(j1=0;j1<2*B;j1++)
        {
            alpha=alpha_j(j1);              
            for(j2=0;j2<2*B;j2++)
            {
                f_j1kj2=func(alpha,beta,gamma_j(j2));
                F_beta[j1][j2]=F_beta[j1][j2]+d_beta_k*(weight[k]*f_j1kj2); 
            }
        }
    }   

    for(j1=0;j1<2*B;j1++)
    {
        F_gamma[j1].init(l_max);
        F_gamma[j1].setZero();
    }

    for(j1=0;j1<2*B;j1++)
    {
        for(j2=0;j2<2*B;j2++)
        {
            gamma=gamma_j(j2);          
            for(l=0;l<=l_max;l++)
            {
                for(m=-l;m<=l;m++)
                {
                    for(n=-l;n<=l;n++)
                    {
                        F_gamma[j1](l,m,n)+=exp(I*((double)n)*gamma)*F_beta[j1][j2](l,m,n);                     
                    }
                }
            }
        }
    }

    F.setZero();
    for(j1=0;j1<2*B;j1++)
    {
        alpha=alpha_j(j1);      
        for(l=0;l<=l_max;l++)
        {
            for(m=-l;m<=l;m++)
            {
                for(n=-l;n<=l;n++)
                {
                    F(l,m,n)+=exp(I*((double)m)*alpha)*F_gamma[j1](l,m,n);                      
                }
            }
        }
    }
   
    return F;
}

fdcl_FFTSO3_matrix_complex fdcl_FFTSO3_complex::forward_transform(std::function <complex<double>(Matrix3)> func) 
{
    // See lambda expression : https://www.geeksforgeeks.org/lambda-expression-in-c/
    return forward_transform([=] (double a, double b, double g)
            {
                return func(Euler3232R(a,b,g));
            });
}

double fdcl_FFTSO3_complex::check_Clebsch_Gordon()
{
    int l1=2, l2=5, l, m1, m2, n1, n2;
    double alpha, beta, gamma;
    fdcl_FFTSO3_matrix_complex D(l1+l2);
    complex<double> y, y_CB={0., 0.};
    double error=0.;

    alpha=(double)rand()/RAND_MAX*2.*M_PI;
    beta=(double)rand()/RAND_MAX*M_PI;
    gamma=(double)rand()/RAND_MAX*2.*M_PI;

    D=wigner_D(alpha,beta,gamma,l1+l2);

    C.compute(l1,l2);

    if(check_verbose)
    {
        cout << "fdcl_FFTSO3_complex::check_Clebsch_Gordon" << endl;
        cout << "l1 = " << l1 << ", l2 = " << l2 << endl;
        cout << "alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << endl;
    }

    for(m1=-l1;m1<=l1;m1++)
        for(n1=-l1;n1<=l1;n1++)
            for(m2=-l2;m2<=l2;m2++)
                for(n2=-l2;n2<=l2;n2++)
                {
                    y=D(l1,m1,n1)*D(l2,m2,n2);
                    y_CB={0., 0.};
                    for(l=max(abs(l1-l2),max(abs(m1+m2),abs(n1+n2))) ;l<=l1+l2; l++)
                        // for(m=-l;m<=l;m++)
                            // for(n=-l;n<=l;n++)
                                y_CB+=C(l,m1+m2,l1,m1,l2,m2)*C(l,n1+n2,l1,n1,l2,n2)*D(l,m1+m2,n1+n2);

                    // cout << "y = " << y << endl;
                    // cout << "y_Clebsch_Gordon = " << y_CB << endl;
                    if(abs(y-y_CB) > error)
                        error = abs(y-y_CB);
                }

    cout << "fdcl_FFTSO3_complex::check_Clebsch_Gordon: error = " << error << endl;

    return error;
}

void fdcl_FFTSO3_complex::check_all()
{
    check_weight();
    check_wigner_d();
    check_deriv_wigner_D();
    check_transform();
    check_Clebsch_Gordon();
    cout << endl;
}

fdcl_FFTSO3_real::fdcl_FFTSO3_real(int l_max)
{
    init(l_max);
}

void fdcl_FFTSO3_real::init(int l_max)
{
    fdcl_FFTSO3_complex::init(l_max);
    U.init(l_max);
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
	T.init(L);
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

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_real::forward_transform(std::function <double(double, double, double)> func)
{
	fdcl_FFTSO3_matrix_complex F_complex(l_max);

	F_complex=fdcl_FFTSO3_complex::forward_transform(func);
	matrix2rsph(l_max);

	for(int l=0; l<=l_max; l++)
		F_complex[l]=T[l]*F_complex[l]*T[l].adjoint();

	return F_complex.real();
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3_real::forward_transform_1(std::function <double(double, double, double)> func)
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
    else if ( (m<0 && n==0) || (m==0 && n < 0) )
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

double fdcl_FFTSO3_real::check_wigner_D_real()
{
    int L=10;
    double alpha, beta, gamma;
    fdcl_FFTSO3_matrix_complex U_2(L);
    fdcl_FFTSO3_matrix_real U(L), U_0(L), U_1(L);
    fdcl_tictoc tictoc;
    double error=1.0;

    alpha=(double)rand()/RAND_MAX*M_PI*2.;
    beta=(double)rand()/RAND_MAX*M_PI;
    gamma=(double)rand()/RAND_MAX*M_PI*2.;

    U=wigner_D_real(alpha,beta,gamma,L);
    U_0=wigner_D_real_0(alpha,beta,gamma,L);
    U_1=wigner_D_real_1(alpha,beta,gamma,L);
    U_2=wigner_D_real_2(alpha,beta,gamma,L);

    if(check_verbose)
    {
        cout << "fdcl_FFTSO3_real::check_wigner_D_real" << endl;
        cout << "alpha=" << alpha << ", beta=" << beta << ", gamma=" << gamma << endl;

        cout << "error from wigner_D_real_0: " << (U-U_0).norm() << endl;
        cout << "error from wigner_D_real_1: " << (U-U_1).norm() << endl;
        cout << "error from wigner_D_real_2: " << (U-U_2.real()).norm() << endl;
    }

    error = (U-U_0).norm() + (U-U_1).norm() + (U-U_2.real()).norm();

    if(check_verbose)
    {
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
    } 

    cout << "fdcl_FFTSO3_real::check_wigner_D_real: error = " << error << endl;
    return error;
}

std::vector<fdcl_FFTSO3_matrix_real> fdcl_FFTSO3_real::deriv_wigner_D_real()
{
    std::vector<fdcl_FFTSO3_matrix_real> u;
    double tmp;
    u.resize(4);
    u[1].init(l_max);
    u[2].init(l_max);
    u[3].init(l_max);   
    
    int l,m,n,p;
    
    for(l=0;l<=l_max;l++)
    {
        for(p=1;p<=l;p++)
        {
            m=p; n=-m;
            u[3](l,m,n)=-m;

            m=-p; n=-m;
            u[3](l,m,n)=-m;

        }
        for(p=2; p<=l; p++)
        {
            m=p; n=m-1;
            tmp=0.5*sqrt((l+m)*(l-m+1));
            u[2](l,m,n)=tmp;

            m=-p; n=m+1;
            u[2](l,m,n)=tmp;

            m=p; n=-m+1;
            u[1](l,m,n)=(m+n)*tmp;

            m=-p; n=-m-1;
            u[1](l,m,n)=(m+n)*tmp;

            m=p-1; n=m+1;
            tmp=0.5*sqrt((l-m)*(l+m+1));
            u[2](l,m,n)=-tmp;

            m=-p+1; n=m-1;
            u[2](l,m,n)=-tmp;

            m=p-1; n=-m-1;
            u[1](l,m,n)=-(m+n)*tmp;

            m=-p+1; n=-m+1;
            u[1](l,m,n)=-(m+n)*tmp;
        }
        if(l>=1)
        {
            tmp=1./sqrt(2.)*sqrt(l*(l+1));
            m=1; n=0;
            u[2](l,m,n)=tmp;

            m=0; n=1;
            u[2](l,m,n)=-tmp;

            m=-1; n=0;
            u[1](l,m,n)=-tmp;

            m=0; n=-1;
            u[1](l,m,n)=tmp;
        }
    }
    
    return u;
}

double fdcl_FFTSO3_real::check_deriv_wigner_D_real()
{
    std::vector<fdcl_FFTSO3_matrix_real> u;
    fdcl_FFTSO3_matrix_real U, U_new;
    std::vector<double>abg;
    Eigen::Matrix<double, 3, 1> ei;
    double eps=1.e-6;
    double error =0.;
    
    u=deriv_wigner_D_real();
    U=wigner_D_real(0,0,0);

    if(check_verbose)
        cout << "fdcl_FFTSO3_real::check_deriv_wigner_D_real" << endl;

    for(int i=1;i<=3;i++)
    {
        if(check_verbose)
            cout << endl << "u_" << i << endl;
        ei.setZero();
        ei(i-1)=1.;
        abg=R2Euler323(expm_SO3(ei*eps));
        U_new=wigner_D_real(abg[0],abg[1],abg[2]);
        for(int l=1;l<=l_max;l++)
        {
            if(check_verbose)
                cout << "l=" << l << ", error= " << ((U_new[l]-U[l])/eps-u[i][l]).norm()/u[i][l].norm() << endl;
            
            error += ((U_new[l]-U[l])/eps-u[i][l]).norm()/u[i][l].norm();
            // cout << u[i][l] << endl << endl; // analytic derivative
            // cout << (U_new[l]-U[l])/eps << endl << endl; // numerical derivative
        }
    }
   
    cout << "fdcl_FFTSO3_real::check_deriv_wigner_D_real: error = " << error*1.e-9 << endl;
    return error*1.e-9;
}

double fdcl_FFTSO3_real::check_Clebsch_Gordon()
{
	int l1=3, l2=5, l, m1, m2, n1, n2;
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


	// c.print();

	// cout << c.c << endl;
    if(check_verbose)
    {
        cout << "fdcl_FFTSO3_real::check_Clebsch_Gordon" << endl;
        cout << "l1 = " << l1 << ", l2 = " << l2 << endl;
        cout << "alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << endl;
    }

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

    
    cout << "fdcl_FFTSO3_real::check_Clebsch_Gordon: error = " << error << endl;

    return error;
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

double fdcl_FFTSO3_real::f_4_check_transform(double alpha, double beta, double gamma)
{
    return inverse_transform(F_4_check,alpha,beta,gamma);
}

double fdcl_FFTSO3_real::check_transform()
{
    double error = 0.;
    F_4_check.init(l_max);
    F_4_check.setRandom();

    auto func= std::bind(&fdcl_FFTSO3_real::f_4_check_transform, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);

    error = (F_4_check-forward_transform(func)).norm();
    cout << "fdcl_FFTSO3_real::check_transform: l_max=" << l_max << " : error = " << error << endl;
    return error;
}

void fdcl_FFTSO3_real::check_all()
{
    check_wigner_D_real();
    check_Clebsch_Gordon();
    check_deriv_wigner_D_real();
    check_transform();
    cout << endl;
}
