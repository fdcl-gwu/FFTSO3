#include <iostream>
#include <vector>
#include <math.h> // pow
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const complex<double> I(0.0,1.0);    

class fdcl_FFTSO3_matrix
{
public:
	int l_max;
	std::vector<MatrixXd> M;
	fdcl_FFTSO3_matrix(){};
	~fdcl_FFTSO3_matrix(){};
	fdcl_FFTSO3_matrix(int l_max); 
	void init(int l_max);
	MatrixXd& operator[](int l); // return l-th matrix 
	double& operator()(int l, int m, int n); // access (m,n)-th element of the l-th matrix
    friend ostream& operator<<(ostream& os, const fdcl_FFTSO3_matrix& M);  	
private:
	void assert_index(int l);
	void assert_index(int l, int m, int n);
};

fdcl_FFTSO3_matrix::fdcl_FFTSO3_matrix(int l_max)
{
	init(l_max);
}

void fdcl_FFTSO3_matrix::init(int l_max)
{
	this->l_max=l_max;
	
	M.resize(l_max+1);
	for(int i=0;i<=l_max;i++)
	{
		M[i].resize(2*i+1,2*i+1);
		M[i].setZero();
	}
}

void fdcl_FFTSO3_matrix::assert_index(int l)
{
	assert(l>=0 && l<=l_max);
}

void fdcl_FFTSO3_matrix::assert_index(int l, int m, int n)
{
	assert_index(l);
	assert(min(m,n) >= -l && max(m,n) <= l);
}

MatrixXd & fdcl_FFTSO3_matrix::operator[](int l)
{
	assert_index(l);
	return M[l];
}

double& fdcl_FFTSO3_matrix::operator()(int l, int m, int n)
{
	assert_index(l,m,n);
	return M[l](m+l,n+l);
}

ostream& operator<<(ostream& os, const fdcl_FFTSO3_matrix& M)  	
{
	for(int l=0;l<=M.l_max;l++)
	{
		os << "l=" << l << endl;
		os << M.M[l] << endl << endl;
	}
	return os;
}




class fdcl_FFTSO3_matrix_complex
{
public:
	int l_max;
	std::vector<MatrixXcd> M;
	fdcl_FFTSO3_matrix_complex(){};
	~fdcl_FFTSO3_matrix_complex(){};
	fdcl_FFTSO3_matrix_complex(int l_max); 
	void init(int l_max);
	MatrixXcd& operator[](int l); // return l-th matrix 
	complex<double>& operator()(int l, int m, int n); // access (m,n)-th element of the l-th matrix
    friend ostream& operator<<(ostream& os, const fdcl_FFTSO3_matrix_complex& M);  	
private:
	void assert_index(int l);
	void assert_index(int l, int m, int n);
};

fdcl_FFTSO3_matrix_complex::fdcl_FFTSO3_matrix_complex(int l_max)
{
	init(l_max);
}

void fdcl_FFTSO3_matrix_complex::init(int l_max)
{
	this->l_max=l_max;
	
	M.resize(l_max+1);
	for(int i=0;i<=l_max;i++)
	{
		M[i].resize(2*i+1,2*i+1);
		M[i].setZero();
	}
}

void fdcl_FFTSO3_matrix_complex::assert_index(int l)
{
	assert(l>=0 && l<=l_max);
}

void fdcl_FFTSO3_matrix_complex::assert_index(int l, int m, int n)
{
	assert_index(l);
	assert(min(m,n) >= -l && max(m,n) <= l);
}

MatrixXcd & fdcl_FFTSO3_matrix_complex::operator[](int l)
{
	assert_index(l);
	return M[l];
}

complex<double>& fdcl_FFTSO3_matrix_complex::operator()(int l, int m, int n)
{
	assert_index(l,m,n);
	return M[l](m+l,n+l);
}

ostream& operator<<(ostream& os, const fdcl_FFTSO3_matrix_complex& M)  	
{
	for(int l=0;l<=M.l_max;l++)
	{
		os << "l=" << l << endl;
		os << M.M[l] << endl << endl;
	}
	return os;
}


class fdcl_FFTSO3
{
public:
	fdcl_FFTSO3_matrix d, D;
	int l_max;
	
	fdcl_FFTSO3(){};
	fdcl_FFTSO3(int);
	~fdcl_FFTSO3(){};
	
	fdcl_FFTSO3_matrix wigner_d(double beta);
	fdcl_FFTSO3_matrix_complex wigner_D(double alpha, double beta, double gamma);
	
	
};

fdcl_FFTSO3::fdcl_FFTSO3(int l_max)
{
	this->l_max=l_max;	
}

fdcl_FFTSO3_matrix fdcl_FFTSO3::wigner_d(double beta)
{
	// M. Blanco and M. Florez and M Bermejo, "Evaluation of the rotation matrices in the basis of real spherical harmonics," Journal of Molecular Structure, 419, pp 19-27, 1997
	fdcl_FFTSO3_matrix d(l_max);
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
	fdcl_FFTSO3_matrix d(l_max);
	fdcl_FFTSO3_matrix_complex D(l_max);
	int l,m,n;
	
	d=wigner_d(beta);

	for(l=0;l<=l_max;l++)
		for(m=-l;m<=l;m++)
			for(n=-l;n<=l;n++)
				D(l,m,n)=d(l,m,n)*exp( -I*(alpha*((double)m) + gamma*((double)n)) );
	
	return D;
}

int main()
{
	
	fdcl_FFTSO3_matrix d(5);
	fdcl_FFTSO3_matrix_complex D(5);
	
	fdcl_FFTSO3 FFTSO3(5);
	
	d=FFTSO3.wigner_d(1.);
	D=FFTSO3.wigner_D(1.,2.,3.);
	cout << D;

	cout << d[3].transpose()*d[3] << endl << endl;
	cout << (D[3].adjoint()*D[3]).real() << endl << endl;

}
