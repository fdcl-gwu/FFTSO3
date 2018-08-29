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

int main()
{
    int l_max=200;
    double dt1, dt2;
    fdcl_tictoc tt;
    fdcl_FFTS2_complex FFTS2(l_max);
    fdcl_FFTS2_real RFFTS2(l_max);
    fdcl_FFTS2_matrix_complex Y(l_max);
    fdcl_FFTSO3_complex FFTSO3(l_max);

    Y.setRandom();

    omp_set_dynamic(0);     // Explicitly disable dynamic teams

    std::vector<double> dt;
    dt.resize(4);
    for(int r=0; r<=2; r++)
    {
        omp_set_num_threads(std::pow(2,r)); // 
        tt.tic();
        for(int i=0; i<=10; i++)
        {
        // FFTS2.nor_assoc_Legendre_poly(0.1,l_max);
        // FFTS2.spherical_harmonics(0.1,0.2,l_max);
        // FFTS2.compute_weight();
        // FFTS2.forward_transform(myfunc, 0);
        // FFTS2.inverse_transform(Y,0.1,0.2);
        // RFFTS2.matrix2rsph(l_max);
        // RFFTS2.spherical_harmonics(0.1,0.2,l_max);
        // RFFTS2.forward_transform(myfunc_real);
        // RFFTS2.inverse_transform(Y.real(),0.1,0.2);
        FFTSO3.wigner_d(0.1,l_max);
        }

        dt[r]=tt.toc();
    }
    cout << "speed up factor " << dt[0]/dt[1] << " " << dt[0]/dt[2] << " " << dt[0]/dt[3] << endl;

    FFTSO3.check_wigner_d();

}
