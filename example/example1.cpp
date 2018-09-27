// elaborated example for real harmonic analysis on SO(3)
#include <iostream>
#include "fdcl_FFTSO3.hpp"

using std::cout;
using std::endl;

int l_max=2;
fdcl::FFTSO3_matrix_real F(l_max); 

// define a real-valued function as a Fourier expansion with random parameter
double func(double alpha, double beta, double gamma)
{
    fdcl::FFTSO3_real RFFTSO3(l_max);
    return RFFTSO3.inverse_transform(F,alpha,beta,gamma);
}

    
int main()
{
    fdcl::FFTSO3_real RFFTSO3(l_max); 

    // Forward transform and inverse transform
    F.setRandom();
    cout << "Check if the forward transform recovers a random Fourier parameter" << endl;
    cout << "error = " << (RFFTSO3.forward_transform(func)-F).norm() << endl << endl;

    // Real harmonics
    double alpha=M_PI/2, beta=M_PI/4, gamma=0;
    cout << "Real harmonics U^l_{m,n} with alpha=" << alpha << ", beta=" << beta << ", gamma=" << gamma << endl;
    cout << "U^l(alpha,beta,gamma)" << RFFTSO3.real_harmonics(alpha,beta,gamma);

    // Derivatives of real harmonics
    cout << "Derivatives of real harmonics" << endl;
    std::vector<fdcl::FFTSO3_matrix_real> u;
    u.resize(4);
    u=RFFTSO3.deriv_real_harmonics();
    cout << "u^l(e_1) = " << u[1] << endl;

    // Clebsch-Gordon coefficients
    int l1=1, l2=2, l=1, m=1, m1=1, m2=-2; 
    cout << "Real Clebsch-Gordon coefficient with (l,m,l1,m1,l2,m2)=(" << l << "," << m << "," << l1 << "," << m1 << "," << l2 << "," << m2 << ")." << endl;
    RFFTSO3.c.compute(l1,l2);
    cout << "c^{l,m}_{l_1,m_1,l_2,m_2} = " << RFFTSO3.c(l,m,l1,m1,l2,m2) << endl;

    return 0;
    
}

