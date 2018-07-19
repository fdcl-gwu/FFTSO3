#include <iostream>
#include <vector>
#include <math.h> // pow
#include <iomanip>      // std::setprecision
#include <Eigen/Dense>
#include <string>

#include "misc_matrix_func.h"
#include "fdcl_FFTSO3_matrix.hpp"

using namespace std;
using namespace Eigen;

typedef Eigen::Matrix<double, 3, 3> Matrix3;

class fdcl_Clebsch_Gordon_matrix
{
    public:
        fdcl_Clebsch_Gordon_matrix(){};
        fdcl_Clebsch_Gordon_matrix(int l1, int l2);
        void init(int l1, int l2);
        ~fdcl_Clebsch_Gordon_matrix(){};
        int l1, l2;
        Eigen::Matrix<double,Dynamic,Dynamic> C;
        Eigen::Matrix<complex<double>,Dynamic,Dynamic> c;

        double& operator()(int l, int m, int l1, int m1, int l2, int m2); 

        int row(int l, int m, int l1, int m1, int l2, int m2);
        int col(int l, int m, int l1, int m1, int l2, int m2);
        void assert_index(int l, int m, int l1, int m1, int l2, int m2);
        void compute_sub(int l, int m, int l1, int l2);
		void compute(int l1, int l2);
        void compute_real(int l1, int l2);
        fdcl_FFTSO3_matrix_complex matrix2rsph(int L);
};

fdcl_Clebsch_Gordon_matrix::fdcl_Clebsch_Gordon_matrix(int l1, int l2)
{
    init(l1,l2);
}

void fdcl_Clebsch_Gordon_matrix::init(int l1, int l2)
{
    int n;
    this->l1=l1;
    this->l2=l2;
    n = (2*l1+1)*(2*l2+1);
    C.resize(n,n);
    C.setZero();
}

int fdcl_Clebsch_Gordon_matrix::row(int l, int m, int l1, int m1, int l2, int m2)
{
    return (l1+m1)*(2*l2+1)+l2+m2;
}

int fdcl_Clebsch_Gordon_matrix::col(int l, int m, int l1, int m1, int l2, int m2)
{
    return l*l-(l2-l1)*(l2-l1)+l+m;
}

void fdcl_Clebsch_Gordon_matrix::assert_index(int l, int m, int l1, int m1, int l2, int m2)
{
    assert( l >= abs(l1-l2) & l <= l1+l2);
    assert( m >= -l & m <=l );
    assert( m1 >= -l1 & m1 <= l1 & m2 >= -l2 & m2 <= l2);
}

double& fdcl_Clebsch_Gordon_matrix::operator() (int l, int m, int l1, int m1, int l2, int m2)
{
    assert_index(l,m,l1,m1,l2,m2);
    return C(row(l,m,l1,m1,l2,m2),col(l,m,l1,m1,l2,m2));
}

void fdcl_Clebsch_Gordon_matrix::compute_sub(int l, int m, int l1, int l2)
{
    // implementation of Straub (2014) Efficient computation of Clebsch-Gordon coefficients http://vixra.org/abs/1403.0263
    int mm, n, x, m1, i;
    double count;
    Eigen::Matrix<double, Dynamic, 1> BB, CC;

    mm = (m-l1-l2+abs(l1-l2+m))/2;
    n = (m+l1+l2-abs(l1-l2-m))/2-mm+1;
    BB.resize(2*n,1);
    CC.resize(n+1,1);
    BB.setZero();
    CC.setZero();

    count=0.;
    CC(n-1)=1.;
    for(x=n-1;x>=1;x--)
    {
        BB(2*x-1) = l1*(l1+1.) + l2*(l2+1.) + 2.*(mm+x)*(m-mm-x) - l*(l+1.);
        BB(2*x-2) = sqrt( (l1*(l1+1.) - (mm+x)*(mm+x-1.))*(l2*(l2+1.) - (m-mm-x)*(m-mm-x+1.)));
        CC(x-1) = -(BB(2*x-1)*CC(x) + BB(2*x)*CC(x+1))/BB(2*x-2);
        count += pow(CC(x-1),2);
    }

    CC(n-1) = sqrt(1./(count+1.));
	for(x=n-1;x>=1;x--)
	{
		CC(x-1) = -(BB(2*x-1)*CC(x) + BB(2*x)*CC(x+1))/BB(2*x-2);
		count += pow(CC(x-1),2);
	}

	i=0;
	for(m1=-l1;m1<=l1;m1++)
	{
		if(abs(m-m1)<=l2)
		{
			C(row(l,m,l1,m1,l2,m-m1),col(l,m,l1,m1,l2,m-m1))=CC(i);
			i+=1;
		}
	}

    // cout << BB << endl << endl;
    // cout << CC << endl;
}

void fdcl_Clebsch_Gordon_matrix::compute(int l1, int l2)
{
    this->init(l1,l2);
	for (int l=abs(l1-l2);l<=l1+l2;l++)
		for (int m=-l;m<=l;m++)
			compute_sub(l,m,l1,l2);

}

fdcl_FFTSO3_matrix_complex fdcl_Clebsch_Gordon_matrix::matrix2rsph(int L)
{
    fdcl_FFTSO3_matrix_complex T(L);
    int l,m;

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

void fdcl_Clebsch_Gordon_matrix::compute_real(int l1, int l2)
{
    fdcl_FFTSO3_matrix_complex T;
    Eigen::Matrix<complex<double>,Dynamic,Dynamic> T12, OTl;
    int n, n1, n2, il;

    n1 = (2*l1+1);
    n2 = (2*l2+1);
    n = n1*n2;
    T.init(l1+l2);
    T = matrix2rsph(l1+l2);
    T12.resize(n,n);
    T12.setZero();
    OTl.resize(n,n);
    OTl.setZero();

    compute(l1,l2);
    for(int i=0;i<2*l1+1;i++)
        for(int j=0;j<2*l1+1;j++)
        {
            T12.block(n2*i,n2*j,n2,n2)=T[l2]*T[l1](i,j);
            // cout << T[l2] << endl;
        }

    il=0;
    for(int l=abs(l1-l2); l<=l1+l2; l++)
    {
        OTl.block(il,il,2*l+1,2*l+1)=T[l];
        il+=2*l+1;
        cout << il << endl;
    }
    c = T12.conjugate()*C*OTl.transpose();

    cout << T[l1] << endl << endl;
    cout << T[l2] << endl << endl;
    cout << T12 << endl << endl;
    cout << T12.conjugate() << endl << endl;
    cout << OTl << endl << endl;
    cout << c*c.adjoint() << endl;
}

int main()
{
    fdcl_Clebsch_Gordon_matrix C;
    int l1, l2, m1, m2, l, m;
    l1=1;
    l2=1;
    m1=1;
    m2=2;

    C.init(l1,l2);

    // for(l=abs(l1-l2);l<=l1+l2;l++)
    // {
        // for(m=-l;m<=l;m++)
        // {
            // cout << C.row(l,m,l1,m1,l2,m2) << ", " << C.col(l,m,l1,m1,l2,m2) << endl;
            // cout << C(l,m,l1,m1,l2,m2) << endl;
        // }
    // }
// 
//
    C.compute(l1,l2);
    C.compute_real(l1,l2);
// cout << C.C.rows() << ", " << C.C.cols() << endl;

    return 0;
}

