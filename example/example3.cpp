#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <string> // for class ETOPO5
#include <fstream>

#include "fdcl_FFTSO3.hpp"
#include "fdcl_FFTS2.hpp"
#include "misc_matrix_func.h"

#include "example3_ETOPO5.hpp"
#include "example3_SSM.hpp"

int main()
{
    int l_max=128;
    fdcl::FFTSO3_real FFTSO3(l_max);
    fdcl::FFTS2_real FFTS2(l_max);
    fdcl::FFTSO3_matrix_real U(l_max);

    fdcl::spherical_shape_matching SSM(l_max);
    fdcl::ETOPO5 ETOPO5;
    std::ofstream fd;
    Eigen::Matrix3d R_true, R_opt;
    bool save_to_file=false;

    //// read ETOPO5 elevation
    ETOPO5.read("./ETOPO5.asc",180,360);

    //// compute Fourier transform of the original ETOPO5 elevation
    auto func= std::bind(&fdcl::ETOPO5::elev, &ETOPO5, std::placeholders::_1, std::placeholders::_2);
    SSM.F = FFTS2.forward_transform([=] (double theta, double phi) {return func(-theta*180./M_PI+90.,phi*180./M_PI);} );
    
    // save the results of Fourier transform for the original elevation
    if(save_to_file)
    {
        fd.open("FFTS2.dat");
        for(int i=0; i<180; i++)
        {
            for(int j=0; j<360; j++)
                fd << FFTS2.inverse_transform(SSM.F,(double)i*M_PI/180., (double)j*M_PI/180.) << " ";

            fd << endl;
        }
        fd.close();
    }
    

    // rotate the original data by R_true, and compute the rotated spherical harmonics
    R_true=fdcl::Euler3232R(M_PI/6., M_PI/3., M_PI/4.);
    U=FFTSO3.real_harmonics(R_true);
    for(int l=0; l<=l_max; l++)
        SSM.G[l]=U[l]*SSM.F[l];


    SSM.G = FFTS2.forward_transform([=] (double theta, double phi) {
            
            Eigen::Vector3d x;
            double theta_rot, phi_rot;
            x << cos(phi)*sin(theta) , sin(phi)*sin(theta), cos(theta);
            x = R_true.transpose()*x;

            theta_rot = acos(x(2));
            if (abs(theta_rot) < 1e-6 || abs(theta_rot-M_PI) < 1e-6)
                phi_rot=0.;
            else
                phi_rot = atan2(x(1)/sin(theta_rot), x(0)/sin(theta_rot));

            if(phi_rot < 0)
                phi_rot += 2*M_PI;

            return func(-theta_rot*180./M_PI+90.,phi_rot*180./M_PI);} );
    
    // save the results of Fourier transform for the rotated elevation
    if(save_to_file)
    {
        fd.open("FFTS2_rot.dat");
        for(int i=0; i<180; i++)
        {
            for(int j=0; j<360; j++)
                fd << FFTS2.inverse_transform(SSM.G,(double)i*M_PI/180., (double)j*M_PI/180.) << " ";

            fd << endl;
        }
        fd.close();
    }

    // perform optimization
    SSM.step_size = 1.e-4;
    R_opt = SSM.opt(fdcl::Euler3232R(0.3,0.3,0.3)); 

    cout << "Ideal cost = " << SSM.J(R_true) << endl;
    cout << "True rotation = " << endl << R_true << endl;
    cout << "Optimized cost = " << SSM.J(R_opt) << endl;
    cout << "Optimized rotation = " << endl << R_opt << endl;

}
