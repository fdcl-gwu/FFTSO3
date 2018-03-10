#include <iostream>
#include <vector>
#include <math.h> // pow
#include <iomanip>      // std::setprecision
#include <Eigen/Dense>
#include "fdcl_FFTSO3_matrix.hpp"
#include "misc_matrix_func.h"

using namespace std;
using namespace Eigen;

typedef Eigen::Matrix<double, 3, 3> Matrix3;


class fdcl_FFTSO3
{
public:
	std::vector<fdcl_FFTSO3_matrix_real> d_beta;
	fdcl_FFTSO3_matrix_complex D, u, F, F0;
	int B, l_max;
	std::vector<double> weight;
	
	fdcl_FFTSO3(){};
	fdcl_FFTSO3(int l_max);
	~fdcl_FFTSO3(){};

	fdcl_FFTSO3_matrix_real wigner_d(double beta, int L);	
	fdcl_FFTSO3_matrix_real wigner_d(double beta);
	fdcl_FFTSO3_matrix_real wigner_d_explicit(double beta);
	
	fdcl_FFTSO3_matrix_complex wigner_D(double alpha, double beta, double gamma);
	fdcl_FFTSO3_matrix_complex wigner_D(double alpha, double beta, double gamma, int L);

	std::vector<fdcl_FFTSO3_matrix_complex> deriv_D();
	complex<double> inverse_transform(fdcl_FFTSO3_matrix_complex, double alpha, double beta, double gamma);
	complex<double> inverse_transform(fdcl_FFTSO3_matrix_complex, Matrix3);
	complex<double> inverse_transform(double alpha, double beta, double gamma);
	complex<double> inverse_transform(Matrix3);
	std::vector<double> character(double beta);
	std::vector<double> compute_weight();
	fdcl_FFTSO3_matrix_complex forward_transform();

	void check_weight();
	void check_wigner_d();
	void check_deriv_D();

	complex<double> f(double alpha, double beta, double gamma);

	Eigen::VectorXd Legendre_poly(double x, int n);
	
private:
	double delta(int ,int );
	double beta_k(int k);
	double alpha_j(int j);
	double gamma_j(int j);
	
	
	
};

fdcl_FFTSO3::fdcl_FFTSO3(int l_max)
{
	this->l_max=l_max;	
	B=l_max+1;
	d_beta.resize(2*B);
	weight.resize(2*B);
	F0.init(l_max);
	F0.setRandom();
}

Eigen::VectorXd fdcl_FFTSO3::Legendre_poly(double x, int N)
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

fdcl_FFTSO3_matrix_real fdcl_FFTSO3::wigner_d_explicit(double beta)
{
	// D. Varshalovich, A. Moskalev, and V. Khersonskii, Quantum Theory of Angular Momentum, World Scientific, 1988, Chapter 4
	
	fdcl_FFTSO3_matrix<double> d(l_max);
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


fdcl_FFTSO3_matrix_real fdcl_FFTSO3::wigner_d(double beta)
{
	return wigner_d(beta,l_max);
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3::wigner_d(double beta, int L)
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
	
	// fill the lower triangualr region
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

fdcl_FFTSO3_matrix_complex fdcl_FFTSO3::wigner_D(double alpha, double beta, double gamma)
{
	return wigner_D(alpha,beta,gamma,l_max);
}

fdcl_FFTSO3_matrix_complex fdcl_FFTSO3::wigner_D(double alpha, double beta, double gamma, int L)
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


std::vector<double> fdcl_FFTSO3::compute_weight()
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


void fdcl_FFTSO3::check_weight()
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
	
	cout << "fdcl_FFTSO3::check_weight" << endl;
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

	
/*	cout << "weight(k)/sin(beta_k)" << endl;
	for(k=0;k<2*B;k++)
		cout << "k=" << k << ": " << weight[k]/sin(beta_k(k)) << endl;
*/		
}

void fdcl_FFTSO3::check_wigner_d()
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
	
	d_beta=wigner_d(beta);

	cout << "fdcl_FFTSO3::check_wigner_d" << endl;
	cout << "matrix orthogonality error: " << (d_beta[l].transpose()*d_beta[l]-I).norm() << endl;
	
	for(int i=0;i<N;i++)
	{
		beta=M_PI/((double)N)*((double)i);
		d_beta=wigner_d(beta);
		d_i_1+=d_beta[l].cwiseProduct(d_beta[l])*sin(beta)*M_PI/((double)N);
		d_i_0+=d_beta[l-1].cwiseProduct(d_beta[l].block(1,1,2*l-1,2*l-1))*sin(beta)*M_PI/((double)N);

	}
	d_i_1*=((double)(2*l+1))/2.;

	cout << "functional orthogonality error: " << endl;
	cout << d_i_0 << endl;
	cout << d_i_1 << endl;

	cout << "difference from the explicit expression" << endl;
	beta=(double)rand()/RAND_MAX*M_PI;
	
	d_beta=wigner_d(beta);
	d_beta_explicit=wigner_d_explicit(beta);
	
	cout << "beta =" << beta << endl;
	cout << "l=0: " << (d_beta[0]-d_beta_explicit[0]).norm() << endl;
	cout << "l=1: " << (d_beta[1]-d_beta_explicit[1]).norm() << endl;
	cout << "l=2: " << (d_beta[2]-d_beta_explicit[2]).norm() << endl;
	cout << "l=3: " << (d_beta[3]-d_beta_explicit[3]).norm() << endl;

	beta=-(double)rand()/RAND_MAX*M_PI;
	
	d_beta=wigner_d(beta);
	d_beta_explicit=wigner_d_explicit(beta);
	
	cout << "beta =" << beta << endl;
	cout << "l=0: " << (d_beta[0]-d_beta_explicit[0]).norm() << endl;
	cout << "l=1: " << (d_beta[1]-d_beta_explicit[1]).norm() << endl;
	cout << "l=2: " << (d_beta[2]-d_beta_explicit[2]).norm() << endl;
	cout << "l=3: " << (d_beta[3]-d_beta_explicit[3]).norm() << endl;

	cout << endl;
	
}

std::vector<double> fdcl_FFTSO3::character(double theta)
{
	std::vector<double> chi;
	fdcl_FFTSO3_matrix_real d(l_max);
		
	chi.resize(l_max+1);
	
	for(int l=0;l<=l_max;l++)
		chi[l]=sin( ((double)2*l+1)/2.*theta )/sin(theta/2.);
	
	return chi;
} 
	
std::vector<fdcl_FFTSO3_matrix_complex> fdcl_FFTSO3::deriv_D()
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

double fdcl_FFTSO3::delta(int i, int j)
{
	double delta=0.;
	
	if (i==j)
		delta=1.;
	
	return delta;
	
}

void fdcl_FFTSO3::check_deriv_D()
{
	std::vector<fdcl_FFTSO3_matrix_complex> u;
	fdcl_FFTSO3_matrix_complex D, D_new;
	std::vector<double>abg;
	Eigen::Matrix<double, 3, 1> ei;
	double eps=1.e-6;
	
	u=deriv_D();
	D=wigner_D(0,0,0);

	cout << "fdcl_FFTSO3::check_deriv_D" << endl;

	for(int i=1;i<=3;i++)
	{
		cout << endl << "u_" << i << endl;
		ei.setZero();
		ei(i-1)=1.;
		abg=R2Euler323(expm_SO3(ei*eps));
		D_new=wigner_D(abg[0],abg[1],abg[2]);

		for(int l=1;l<=l_max;l++)
		{
			cout << "l=" << l << ", error= " << ((D_new[l]-D[l])/eps-u[i][l]).norm()/u[i][l].norm() << endl;
//			cout << u[i][l] << endl << endl; // analytic derivative
//			cout << (D_new[l]-D[l])/eps << endl; // numerical derivative
		}
	}

}

complex<double> fdcl_FFTSO3::inverse_transform(fdcl_FFTSO3_matrix_complex F, double alpha, double beta, double gamma)
{
	complex<double> f=0.;
	int l,m,n;
	fdcl_FFTSO3_matrix_complex D;
	
	D=wigner_D(alpha,beta,gamma);	
	
	for(l=0;l<=l_max;l++)
		for(m=-l;m<=l;m++)
			for(n=-l;n<=l;n++)
				f+=((double) 2*l+1 )* F(l,m,n)*D(l,m,n);
	
	return f;
}

complex<double> fdcl_FFTSO3::inverse_transform(double alpha, double beta, double gamma)
{
	return inverse_transform(this->F, alpha, beta, gamma);
}

complex<double> fdcl_FFTSO3::inverse_transform(fdcl_FFTSO3_matrix_complex F, Matrix3 R)
{
	std::vector<double> abg;
	
	abg.resize(3);
	abg=R2Euler323(R);
	
	return inverse_transform(F,abg[0],abg[1],abg[2]);
}

complex<double> fdcl_FFTSO3::inverse_transform(Matrix3 R)
{
	std::vector<double> abg;
	
	abg.resize(3);
	abg=R2Euler323(R);
	
	return inverse_transform(abg[0],abg[1],abg[2]);	
}

double fdcl_FFTSO3::beta_k(int k)
{
	return ((double)(2*k+1))*M_PI/4./((double)B);
}

double fdcl_FFTSO3::alpha_j(int j)
{
	return ((double)j)*M_PI/((double)B);
}

double fdcl_FFTSO3::gamma_j(int j)
{
	return alpha_j(j);
}

complex<double> fdcl_FFTSO3::f(double alpha, double beta, double gamma)
{
/*	fdcl_FFTSO3_matrix_complex D(l_max);
	complex<double> y={0.,0.};
	
	D=wigner_D(alpha,beta,gamma);
	
	for(int l=0;l<=l_max;l++)
		for(int m=-l;m<=l;m++)
			for(int n=-l;n<=l;n++)
				y+=((double)2*l+1)*D(l,m,n);
	
	
	return y;
*/
		
//	return Euler3232R(alpha,beta,gamma).trace();
	
	return alpha+beta+gamma;
//	return inverse_transform(F0,alpha,beta,gamma);
	
}



fdcl_FFTSO3_matrix_complex fdcl_FFTSO3::forward_transform()
{
	fdcl_FFTSO3_matrix_complex F_beta[2*l_max][2*l_max];
	int j1, j2, k, l, m, n;
	fdcl_FFTSO3_matrix_real d_beta_k(l_max);
	fdcl_FFTSO3_matrix_complex F_gamma[2*l_max];
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
				f_j1kj2=f(alpha,beta,gamma);
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
		
/*	for(j1=0;j1<2*l_max;j1++)
	{
		for(j2=0;j2<2*l_max;j2++)
		{
			F_beta[j1][j2].init(l_max);
			F_beta[j1][j2].setZero();
		}
	}

	compute_weight();
	
	for(k=0;k<2*l_max;k++)
	{
		d_beta_k=wigner_d(beta_k(k));
		for(l=0;l<=l_max;l++)
		{
			for(j1=0;j1<2*l_max;j1++)
			{
				for(j2=0;j2<2*l_max;j2++)
				{
					F_beta[j1][j2][l]+=weight[k]*f(alpha_j(j1),beta_k(k),gamma_j(j2))*d_beta_k[l];
				}
			}
		}
		
			
	}
	
	
	for(j1=0;j1<2*l_max;j1++)
	{
		F_gamma[j1].init(l_max);
		F_gamma[j1].setZero();
	}

	for(j1=0;j1<2*l_max;j1++)
	{
		for(j2=0;j2<2*l_max;j2++)
		{
			for(l=0;l<=l_max;l++)
			{
				for(m=-l;m<=l;m++)
				{
					for(n=-l;n<=l;n++)
					{
						F_gamma[j1](l,m,n)+=exp(I*((double)n)*gamma_j(j2))*F_beta[j1][j2](l,m,n);						
					}
				}
			}
		}
	}
	

	F.setZero();
	for(j1=0;j1<2*l_max;j1++)
	{
		for(l=0;l<=l_max;l++)
		{
			for(m=-l;m<=l;m++)
			{
				for(n=-l;n<=l;n++)
				{
					F(l,m,n)+=exp(I*((double)m)*alpha_j(j1))*F_gamma[j1](l,m,n);						
				}
			}
		}
	}
*/	
	return F;
}

int main()
{
	int l_max=7;
	fdcl_FFTSO3_matrix_real d(l_max), d1(l_max);
	fdcl_FFTSO3_matrix_complex D(l_max), F(l_max);

	fdcl_FFTSO3 FFTSO3(l_max);

	FFTSO3.check_weight();
//	FFTSO3.check_wigner_d();

	F=FFTSO3.forward_transform();

//	for(int l=0;l<=l_max;l++)
//		cout << (FFTSO3.F0[l]-F[l]).norm() << endl;
	
	cout << FFTSO3.f(1.,2.,3.) << endl;
	cout << FFTSO3.inverse_transform(F,1.,2.,3.) << endl;
	

}
