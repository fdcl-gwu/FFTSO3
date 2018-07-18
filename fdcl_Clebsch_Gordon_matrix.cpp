#include "fdcl_Clebsch_Gordon_matrix.hpp"

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

