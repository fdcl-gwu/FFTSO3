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

complex<double> myfuncR(double alpha, double beta, double gamma)
{
    return cos(alpha)*sin(beta)*sin(gamma);
}
int main()
{
    int l_max=50;
    double dt1, dt2;
    fdcl_tictoc tt;
    fdcl_FFTS2_complex FFTS2(l_max);
    fdcl_FFTS2_real RFFTS2(l_max);
    fdcl_FFTS2_matrix_complex Y(l_max);
    fdcl_FFTSO3_complex FFTSO3(l_max);
    fdcl_FFTSO3_matrix_complex F(l_max);

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
            FFTSO3.inverse_transform(F,0.1,0.2,0.3);

        }

        dt[r]=tt.toc();
    }
    cout << "speed up factor " << dt[0]/dt[1] << " " << dt[0]/dt[2] << " " << dt[0]/dt[3] << endl;

    FFTSO3.check_wigner_d();
    FFTSO3.init(5);
    FFTSO3.check_transform();

    FFTS2.init(4);
    FFTS2.check_transform();

    FFTS2.init(3);
    FFTS2.check_transform();


}
