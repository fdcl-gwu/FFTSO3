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

fdcl_Clebsch_Gordon_real::fdcl_Clebsch_Gordon_real(int l1, int l2)
{
    init(l1,l2);
}

void fdcl_Clebsch_Gordon_matrix::compute(int l1, int l2)
{
    init(l1,l2);
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

void fdcl_Clebsch_Gordon_real::init(int l1, int l2)
{
    fdcl_Clebsch_Gordon_matrix::init(l1,l2);
    int n = (2*l1+1)*(2*l2+1);
    c.resize(n,n);
    c.setZero();
}


void fdcl_Clebsch_Gordon_real::print()
{
    for(int m1=-l1; m1<=l1; m1++)
        for(int m2=-l2; m2<=l2; m2++)
            for(int l=abs(l1-l2);l<=l1+l2;l++)
                for(int m=-l; m<=l; m++)
                {
                    if(abs(c(row(l,m,l1,m1,l2,m2),col(l,m,l1,m1,l2,m2))) > 1e-5)
                    {
                        cout << "l,m,m1,m2 " << l << " " << m << " " << m1 << " " << m2 <<   endl;
                        cout << c(row(l,m,l1,m1,l2,m2),col(l,m,l1,m1,l2,m2)) << endl;
                    }
                }
}

void fdcl_Clebsch_Gordon_real::compute_0(int l1, int l2)
{
    // matrix form: slower
    init(l1,l2);
    fdcl_FFTSO3_matrix_complex T;
    T.init(l1+l2);
    T = matrix2rsph(l1+l2);
    std::vector<int> P1, P2;

    fdcl_Clebsch_Gordon_matrix::compute(l1,l2);

    for(int m1=-l1; m1<=l1; m1++)
        for(int m2=-l2; m2<=l2; m2++)
            for(int l=abs(l1-l2); l<=l1+l2; l++)
                for(int m=-l; m<=l; m++)
                {
                    P1.clear();
                    if(m1==0)
                        P1.insert(P1.end(),{0});
                    else
                        P1.insert(P1.end(),{-m1,m1});


                    P2.clear();
                    if(m2==0)
                        P2.insert(P2.end(),{0});
                    else
                        P2.insert(P2.end(),{-m2,m2});

                    for(int p1 : P1)
                        for(int p2 : P2)
                            if(abs(p1+p2) == abs(m))
                                c(row(l,m,l1,m1,l2,m2),col(l,m,l1,m1,l2,m2))+=conj(T(l1,m1,p1))*conj(T(l2,m2,p2))*C(row(l,p1+p2,l1,p1,l2,p2),col(l,p1+p2,l1,p1,l2,p2))*T(l,m,p1+p2);

                }
    // alternative method: matrix computation using Kronecker product : slower
    // int N, N1, N2;
    // int il;
    // Eigen::Matrix<complex<double>,Dynamic,Dynamic> T12, OTl, c_new;
    // N1 = (2*l1+1);
    // N2 = (2*l2+1);
    // N = N1*N2;
    // T12.resize(N,N);
    // T12.setZero();
    // OTl.resize(N,N);
    // OTl.setZero();
    // c_new.resize(N,N);
    // c_new.setZero();
// 
    // for(int i=0;i<2*l1+1;i++)
        // for(int j=0;j<2*l1+1;j++)
            // T12.block(N2*i,N2*j,N2,N2)=T[l2]*T[l1](i,j);
// 
    // il=0;
    // for(int l=abs(l1-l2); l<=l1+l2; l++)
    // {
        // OTl.block(il,il,2*l+1,2*l+1)=T[l];
        // il+=2*l+1;
    // }
    // c_new = T12.conjugate()*C*OTl.transpose();
// 
    // cout << "c error = " << (c-c_new).norm() << endl;
    

    // alternative method: element-wise computation with triple summation: slowest
    // for(int m1=-l1; m1<=l1; m1++)
        // for(int m2=-l2; m2<=l2; m2++)
            // for(int l=abs(l1-l2); l<=l1+l2; l++)
                // for(int m=-l; m<=l; m++)
                // {
                    // for(p1=-l1; p1<=l1; p1++)
                        // for(p2=-l2; p2<=l2; p2++)
                            // for(int p=-l; p<=l; p++)
                                // c_new(row(l,m,l1,m1,l2,m2),col(l,m,l1,m1,l2,m2))+=conj(T(l1,m1,p1))*conj(T(l2,m2,p2))*C(row(l,p,l1,p1,l2,p2),col(l,p,l1,p1,l2,p2))*T(l,m,p);
                // }}
}


void fdcl_Clebsch_Gordon_real::compute(int l1, int l2)
{
    int m1, m2, p1, p2;
    double one_over_sqrt8=1./sqrt(8.);

    init(l1,l2);
    fdcl_Clebsch_Gordon_matrix::compute(l1,l2);


    // case 0
    for(int l=abs(l1-l2); l<=l1+l2; l++)
        c(row(l,0,l1,0,l2,0),col(l,0,l1,0,l2,0))=C(row(l,0,l1,0,l2,0),col(l,0,l1,0,l2,0));

    // case 1,9
    for (p1=1; p1<=l1; p1++)
    {
        m1=p1;
        compute_sub_01(l1,m1,l2,0,{0.5,0.},{0.,0.5});  
        m1=-p1;
        compute_sub_01(l1,m1,l2,0,{0.5,0.},{0.,0.5});  
    }

    // case 5,13
    for (p2=1; p2<=l2; p2++)
    {
        m2=p2;
        compute_sub_01(l1,0,l2,m2,{0.5,0.},{0.,0.5});  
        m2=-p2;
        compute_sub_01(l1,0,l2,m2,{0.5,0.},{0.,0.5});  
    }

    // case 3, 7, 11, 15
    for (p1=1; p1<=min(l1,l2); p1++)
    {
        // case 3
        m1=p1; m2=m1;
        compute_sub_01(l1,m1,l2,m2,{one_over_sqrt8,0.}, {0., one_over_sqrt8});  
        compute_sub_23(l1,m1,l2,m2,{pow(-1.,m2)/2.,0}, {pow(-1.,m2)/2.,0});  

        // case 7
        m1=-p1; m2=p1;
        compute_sub_01(l1,m1,l2,m2,{0., pow(-1.,m1)/2.}, {0., pow(-1.,m1)/2.});
        compute_sub_23(l1,m1,l2,m2,{one_over_sqrt8,0}, {0., one_over_sqrt8});

        // case 11
        m1=-p1; m2=-p1;
        compute_sub_01(l1,m1,l2,m2,{0., one_over_sqrt8}, {-one_over_sqrt8, 0.});
        compute_sub_23(l1,m1,l2,m2,{pow(-1.,m1)/2.,0}, {pow(-1.,m1)/2., 0.});

        // case 15
        m1=p1; m2=-p1;
        compute_sub_01(l1,m1,l2,m2,{0., pow(-1.,m1)/2.}, {0., pow(-1.,m1)/2.});  
        compute_sub_23(l1,m1,l2,m2,{0., -one_over_sqrt8}, {one_over_sqrt8, 0.});  
    }

    // case 2, 8, 10, 16
    for (p1=2; p1<=l1; p1++)
        for (p2=1; p2 <= min(m1-1, l2); p2++)
        {    
            // case 2
            m1=p1; m2=p2;
            compute_sub_01(l1,m1,l2,m2,{one_over_sqrt8,0.}, {0., one_over_sqrt8});  
            compute_sub_23(l1,m1,l2,m2,{pow(-1.,m2)*one_over_sqrt8,0}, {0, pow(-1.,m2)*one_over_sqrt8});  

            // case 8
            m1=-p1; m2=p2;
            compute_sub_01(l1,m1,l2,m2,{pow(-1.,m2)*one_over_sqrt8,0.}, {0., pow(-1.,m2)*one_over_sqrt8});  
            compute_sub_23(l1,m1,l2,m2,{one_over_sqrt8,0}, {0,one_over_sqrt8});  

            // case 10
            m1=-p1; m2=-p2;
            compute_sub_01(l1,m1,l2,m2,{0., one_over_sqrt8}, {-one_over_sqrt8, 0.});  
            compute_sub_23(l1,m1,l2,m2,{0., -pow(-1.,m2)*one_over_sqrt8}, {pow(-1.,m2)*one_over_sqrt8, 0.});  

            // case 16
            m1=p1; m2=-p2;
            compute_sub_01(l1,m1,l2,m2,{0., pow(-1.,m2)*one_over_sqrt8}, {-pow(-1.,m2)*one_over_sqrt8, 0.});  
            compute_sub_23(l1,m1,l2,m2,{0., -one_over_sqrt8}, {one_over_sqrt8, 0.});  
        }

    // case 4, 6, 12, 14
    for (p1=1; p1<=min(l1, l2-1); p1++)
        for (p2=p1+1; p2<=l2; p2++)
        {    
            // case 4
            m1=p1; m2=p2;
            compute_sub_01(l1,m1,l2,m2,{one_over_sqrt8,0.}, {0., one_over_sqrt8});  
            compute_sub_23(l1,m1,l2,m2,{0., -pow(-1.,m1)*one_over_sqrt8}, {pow(-1.,m1)*one_over_sqrt8, 0.});  

            // case 6
            m1=-p1; m2=p2;
            compute_sub_01(l1,m1,l2,m2,{0., pow(-1.,m1)*one_over_sqrt8}, {-pow(-1.,m1)*one_over_sqrt8, 0.});  
            compute_sub_23(l1,m1,l2,m2,{one_over_sqrt8,0}, {0, one_over_sqrt8});  

            // case 12
            m1=-p1; m2=-p2;
            compute_sub_01(l1,m1,l2,m2,{0., one_over_sqrt8}, {-one_over_sqrt8, 0.});
            compute_sub_23(l1,m1,l2,m2,{pow(-1.,m1)*one_over_sqrt8,0.}, {0., pow(-1.,m1)*one_over_sqrt8});  

            // case 14
            m1=p1; m2=-p2;
            compute_sub_01(l1,m1,l2,m2,{pow(-1.,m1)*one_over_sqrt8,0.}, {0., pow(-1.,m1)*one_over_sqrt8});  
            compute_sub_23(l1,m1,l2,m2,{0., -one_over_sqrt8}, {one_over_sqrt8, 0.});  
        }

}

complex<double>& fdcl_Clebsch_Gordon_real::operator() (int l, int m, int l1, int m1, int l2, int m2)
{
    assert_index(l,m,l1,m1,l2,m2);
    return c(row(l,m,l1,m1,l2,m2),col(l,m,l1,m1,l2,m2));
}

void fdcl_Clebsch_Gordon_real::compute_sub_01(int l1, int m1, int l2, int m2, complex<double> ratio0, complex<double> ratio1)
{
    double eta, zeta, C_complex;
    int l,m;

    for(l=max(abs(l1-l2),abs(m1+m2)); l<=l1+l2; l++)
    {
        eta=pow(-1.,l1+l2-l)+1.;
        zeta=pow(-1.,l1+l2-l)-1.;
        m=m1+m2;
        C_complex=C(row(l,m,l1,m1,l2,m2),col(l,m,l1,m1,l2,m2));

        c(row(l,m,l1,m1,l2,m2),col(l,m,l1,m1,l2,m2))=(std::real(ratio0)*eta+std::imag(ratio0)*I*zeta)*C_complex;

        m=-m;
        c(row(l,m,l1,m1,l2,m2),col(l,m,l1,m1,l2,m2))=(std::real(ratio1)*eta+std::imag(ratio1)*I*zeta)*C_complex;
    }
}

void fdcl_Clebsch_Gordon_real::compute_sub_23(int l1, int m1, int l2, int m2, complex<double>ratio2, complex<double>ratio3)
{
    double eta, zeta, C_complex;
    int l,m;

    for(l=max(abs(l1-l2),abs(m1-m2)); l<=l1+l2; l++)
    {
        eta=pow(-1.,l1+l2-l)+1.;
        zeta=pow(-1.,l1+l2-l)-1.;
        m=m1-m2;
        C_complex=C(row(l,m,l1,m1,l2,-m2),col(l,m,l1,m1,l2,-m2));

        c(row(l,m,l1,m1,l2,m2),col(l,m,l1,m1,l2,m2))=(std::real(ratio2)*eta+std::imag(ratio2)*I*zeta)*C_complex;
        m=-m;
        c(row(l,m,l1,m1,l2,m2),col(l,m,l1,m1,l2,m2))=(std::real(ratio3)*eta+std::imag(ratio3)*I*zeta)*C_complex;
    }
}
