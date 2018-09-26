#include <iostream>
#include "fdcl_FFTSO3.hpp"

// define a real-valued function on SO(3)
double func(Eigen::Matrix3d R)
{
    return R.trace();
}
    
int main()
{
    int l_max=2;  // the maximum order of Fourier transform
    fdcl::FFTSO3_real RFFTSO3(l_max); // FFTSO3_real object for real harmonic analysis on SO(3)
    fdcl::FFTSO3_matrix_real F(l_max); // FFTSO3_matrix_real object to save real valued Fourier parameters
    Eigen::Matrix3d R; 

    F=RFFTSO3.forward_transform(func); // perform forward transform
    std::cout << "Fourier parameter" << std::endl << F << std::endl; // show Fourier parameters

    R.setIdentity(); // R is set to the identity matrix
    std::cout << RFFTSO3.inverse_transform(F,R) << std::endl; // compute the inverse transform at the identity and print it

    return 0;
}

