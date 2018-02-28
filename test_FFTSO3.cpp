#include <iostream>
#include <vector>
#include <math.h> // pow
#include <iomanip>      // std::setprecision
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const complex<double> I(0.0,1.0);    

template <class ScalarType>
class fdcl_FFTSO3_matrix
{
public:
	int l_max;
	std::vector<Eigen::Matrix<ScalarType,Dynamic,Dynamic>> M;
	fdcl_FFTSO3_matrix(){};
	~fdcl_FFTSO3_matrix(){};
	fdcl_FFTSO3_matrix(int l_max); 
	void init(int l_max);
	Eigen::Matrix<ScalarType,Dynamic,Dynamic>& operator[](int l); // return l-th matrix 
	ScalarType& operator()(int l, int m, int n); // access (m,n)-th element of the l-th matrix
	
	template<typename _ScalarType>
    friend ostream& operator<<(ostream& os, const fdcl_FFTSO3_matrix<_ScalarType>& M);  	
private:
	void assert_index(int l);
	void assert_index(int l, int m, int n);
};

template <class ScalarType>
fdcl_FFTSO3_matrix<ScalarType>::fdcl_FFTSO3_matrix(int l_max)
{
	init(l_max);
}

template <class ScalarType>
void fdcl_FFTSO3_matrix<ScalarType>::init(int l_max)
{
	this->l_max=l_max;
	
	M.resize(l_max+1);
	for(int i=0;i<=l_max;i++)
	{
		M[i].resize(2*i+1,2*i+1);
		M[i].setZero();
	}
}

template <class ScalarType>
void fdcl_FFTSO3_matrix<ScalarType>::assert_index(int l)
{
	assert(l>=0 && l<=l_max);
}

template <class ScalarType>
void fdcl_FFTSO3_matrix<ScalarType>::assert_index(int l, int m, int n)
{
	assert_index(l);
	assert(min(m,n) >= -l && max(m,n) <= l);
}

template <class ScalarType>
Matrix<ScalarType,Dynamic,Dynamic> & fdcl_FFTSO3_matrix<ScalarType>::operator[](int l)
{
	assert_index(l);
	return M[l];
}

template <class ScalarType>
ScalarType& fdcl_FFTSO3_matrix<ScalarType>::operator()(int l, int m, int n)
{
	assert_index(l,m,n);
	return M[l](m+l,n+l);
}

template <class ScalarType>
ostream& operator<< (ostream& os, const fdcl_FFTSO3_matrix<ScalarType>& M)
{
	for(int l=0;l<=M.l_max;l++)
	{
		os << "l=" << l << endl;
		os << M.M[l] << endl << endl;
	}
	return os;
}


typedef fdcl_FFTSO3_matrix<double> fdcl_FFTSO3_matrix_real;
typedef fdcl_FFTSO3_matrix<complex<double>> fdcl_FFTSO3_matrix_complex;

class fdcl_FFTSO3
{
public:
	std::vector<fdcl_FFTSO3_matrix_real> d_beta;
	fdcl_FFTSO3_matrix_complex D;
	int l_max;
	std::vector<double> weight;
	
	fdcl_FFTSO3(){};
	fdcl_FFTSO3(int l_max);
	~fdcl_FFTSO3(){};
	
	fdcl_FFTSO3_matrix_real wigner_d(double beta);
	fdcl_FFTSO3_matrix_complex wigner_D(double alpha, double beta, double gamma);

	void compute_d_beta();
	void compute_weight();
	double beta_k(int k);
	void check_weight();
private:

	
	
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
int main()
{
	int l_max=20;
	fdcl_FFTSO3_matrix_real d(l_max);
	fdcl_FFTSO3_matrix_complex D(l_max);

	fdcl_FFTSO3 FFTSO3(l_max);
	d=FFTSO3.wigner_d(1.5);
	D=FFTSO3.wigner_D(1.5,2.3,-5.);
	
	cout << d[3].transpose()*d[3] << endl << endl;
	cout << (D[10].adjoint()*D[10]).real() << endl << endl;
	
	FFTSO3.compute_weight();
	
	cout << setprecision(12);
	for(int k=0;k<2*l_max;k++)
		cout << FFTSO3.weight[k] << endl;
	
	FFTSO3.check_weight();
}
