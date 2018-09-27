// elaborated example for complex harmonic analysis on SO(3)
#include <iostream>
#include "fdcl_FFTSO3.hpp"

using std::cout;
using std::endl;

int l_max=2;
fdcl::FFTSO3_matrix_complex F(l_max); 

// define a complex-valued function as a Fourier expansion with random parameter
std::complex<double> func(double alpha, double beta, double gamma)
{
    fdcl::FFTSO3_complex CFFTSO3(l_max);
    return CFFTSO3.inverse_transform(F,alpha,beta,gamma);
}

    
int main()
{
    fdcl::FFTSO3_complex CFFTSO3(l_max); 

    // Forward transform and inverse transform
    F.setRandom();
    cout << "Check if the forward transform recovers a random Fourier parameter" << endl;
    cout << "error = " << (CFFTSO3.forward_transform(func)-F).norm() << endl << endl;

    // Real harmonics
    double alpha=M_PI/2, beta=M_PI/4, gamma=0;
    cout << "Wigner D-matrix with alpha=" << alpha << ", beta=" << beta << ", gamma=" << gamma << endl;
    cout << "D^l(alpha,beta,gamma) = " << endl << CFFTSO3.wigner_D(alpha,beta,gamma);

    // Derivatives of real harmonics
    cout << "Derivatives of Wigner D-matrix" << endl;
    std::vector<fdcl::FFTSO3_matrix_complex> u;
    u.resize(4);
    u=CFFTSO3.deriv_wigner_D();
    cout << "u^l(e_1) = " << endl << u[1] << endl;

    // Clebsch-Gordon coefficients
    int l1=1, l2=2, l=1, m=1, m1=1, m2=-2; 
    cout << "Clebsch-Gordon coefficient with (l,m,l1,m1,l2,m2)=(" << l << "," << m << "," << l1 << "," << m1 << "," << l2 << "," << m2 << ")." << endl;
    CFFTSO3.C.compute(l1,l2);
    cout << "C^{l,m}_{l_1,m_1,l_2,m_2} =" << CFFTSO3.C(l,m,l1,m1,l2,m2) << endl;

    return 0;
    
}

