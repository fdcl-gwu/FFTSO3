#include <iostream>
#include <vector>
#include <math.h> // pow
#include <iomanip>      // std::setprecision
#include <Eigen/Dense>
#include "fdcl_FFTSO3_matrix.hpp"
#include "misc_matrix_func.h"

using namespace std;
using namespace Eigen;

class fdcl_FFTSO3
{
public:
	std::vector<fdcl_FFTSO3_matrix_real> d_beta;
	fdcl_FFTSO3_matrix_complex D, u;
	int l_max;
	std::vector<double> weight;
	
	fdcl_FFTSO3(){};
	fdcl_FFTSO3(int l_max);
	~fdcl_FFTSO3(){};
	
	fdcl_FFTSO3_matrix_real wigner_d(double beta);
	fdcl_FFTSO3_matrix_complex wigner_D(double alpha, double beta, double gamma);
	std::vector<fdcl_FFTSO3_matrix_complex> deriv_D();

	void compute_d_beta();
	void compute_weight();
	double beta_k(int k);
	void check_weight();
	void check_wigner_d();
	std::vector<double> character(double beta);
	void check_deriv_D();
private:
	double delta(int ,int );
	
	
};

fdcl_FFTSO3::fdcl_FFTSO3(int l_max)
{
	this->l_max=l_max;	
	d_beta.resize(2*l_max);
	weight.resize(2*l_max);
}

fdcl_FFTSO3_matrix_real fdcl_FFTSO3::wigner_d(double beta)
{
	// M. Blanco and M. Florez and M Bermejo, "Evaluation of the rotation matrices in the basis of real spherical harmonics," Journal of Molecular Structure, 419, pp 19-27, 1997
	fdcl_FFTSO3_matrix<double> d(l_max);
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
	for(l=2;l<=l_max;l++)
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
	for(l=1;l<=l_max;l++)
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
	fdcl_FFTSO3_matrix_real d(l_max);
	fdcl_FFTSO3_matrix_complex D(l_max);
	int l,m,n;
	
	d=wigner_d(beta);

	for(l=0;l<=l_max;l++)
		for(m=-l;m<=l;m++)
			for(n=-l;n<=l;n++)
				D(l,m,n)=d(l,m,n)*exp( -I*(alpha*((double)m) + gamma*((double)n)) );
	
	return D;
}

void fdcl_FFTSO3::compute_d_beta()
{
	for(int k=0;k<2*l_max;k++)
		d_beta[k]=wigner_d(beta_k(k));
}

void fdcl_FFTSO3::compute_weight()
{	
    int j, k;
    double factor;
    double sum;

    factor = M_PI/((double)(4*l_max)) ;

	for(j=0;j<2*l_max;j++)
	{
		sum=0.0;
        for(k=0;k<l_max;k++)
			sum+=1./((double)(2*k+1))*sin((double)((2*j+1)*(2*k+1))*factor);
		
        sum*=1./((double)l_max)*sin((double)(2*j+1)*factor);
      
        weight[j]=sum;
	}
}

double fdcl_FFTSO3::beta_k(int k)
{
	return ((double)(2*k+1))*M_PI/4./((double)l_max);
}

void fdcl_FFTSO3::check_weight()
{
	fdcl_FFTSO3_matrix<double> d(l_max);
	std::vector<double> sum;
	sum.resize(l_max+1);
	for (int l=0;l<=l_max;l++)
		sum[l]=0.;


	this->compute_weight();
	
	for (int k=0;k<2*l_max; k++)
	{
		d=wigner_d(beta_k(k));
		for(int l=0;l<=l_max;l++)
		{
			sum[l]+=d(l,0,0)*weight[k];
		}
	}
	
	cout << "fdcl_FFTSO3::check_weight" << endl;
	for (int l=0;l<=l_max;l++)
		cout << sum[l] << endl;
		
}

void fdcl_FFTSO3::check_wigner_d()
{
	int l=3, N=1000;
	double beta;
	MatrixXd I, d_i_0, d_i_1;
	fdcl_FFTSO3_matrix_real d_beta(l);
	
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
int main()
{
	int l_max=100;
	fdcl_FFTSO3_matrix_real d(l_max), d1(l_max);
	fdcl_FFTSO3_matrix_complex D(l_max);

	std::vector<fdcl_FFTSO3_matrix_complex> u;

//	fdcl_FFTSO3 FFTSO3(l_max);
	
//	FFTSO3.check_deriv_D();
}
