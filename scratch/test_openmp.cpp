#include <omp.h>
#include <iostream>
#include <vector>

#include "fdcl_tictoc.hpp"
#include "fdcl_FFTS2.hpp"
#include "fdcl_FFTSO3.hpp"

using namespace std;

double myfunc_real(double theta, double phi)
{
    return cos(theta)*sin(phi);
}
complex<double> myfunc(double theta, double phi)
{
    return cos(theta)*sin(phi);
}

double myfuncR_real(double alpha, double beta, double gamma)
{
    return cos(alpha)*sin(beta)*sin(gamma);
}
complex<double> myfuncR(double alpha, double beta, double gamma)
{
    return cos(alpha)*sin(beta)*sin(gamma);
}
int main()
{
    int l_max=50;
    fdcl::tictoc tt;
    fdcl::FFTS2_complex FFTS2(l_max);
    fdcl::FFTS2_real RFFTS2(l_max);
    fdcl::FFTS2_matrix_complex Y(l_max);
    fdcl::FFTSO3_complex FFTSO3(l_max);
    fdcl::FFTSO3_real RFFTSO3(l_max);
    fdcl::FFTSO3_matrix_complex F(l_max);
    fdcl::FFTSO3_matrix_real F_real(l_max);

    Y.setRandom();
    F.setRandom();

    omp_set_dynamic(0);     // Explicitly disable dynamic teams

    // omp_set_num_threads(1);
    // cout << FFTSO3.forward_transform(myfuncR) << endl;
// 
    // omp_set_num_threads(2);
    // cout << FFTSO3.forward_transform(myfuncR) << endl;


    std::vector<double> dt;
    dt.resize(4);
    for(int r=0; r<=3; r++)
    {
        omp_set_num_threads(std::pow(2,r)); // 
        tt.tic();
        for(int i=0; i<=0; i++)
        {
            // FFTSO3.forward_transform(myfuncR);
            // FFTSO3.inverse_transform(F,0.1,0.2,0.3);
            // RFFTSO3.forward_transform(myfuncR_real);
            // RFFTSO3.real_harmonics(0.1,0.2,0.3,l_max);
            // RFFTSO3.inverse_transform(F_real,0.1,0.2,0.3);
            FFTS2.forward_transform(myfunc);

        }

        dt[r]=tt.toc();
    }
    cout << "speed up factor " << dt[0]/dt[1] << " " << dt[0]/dt[2] << " " << dt[0]/dt[3] << endl;

    RFFTS2.init(4);
    RFFTS2.check_transform();
    RFFTS2.init(5);
    RFFTS2.check_transform();

    FFTS2.init(4);
    FFTS2.check_transform();
    FFTS2.init(5);
    FFTS2.check_transform();

    FFTSO3.init(4);
    FFTSO3.check_transform();
    FFTSO3.init(5);
    FFTSO3.check_transform();

    RFFTSO3.init(4);
    RFFTSO3.check_transform();
    RFFTSO3.init(5);
    RFFTSO3.check_transform();
}
