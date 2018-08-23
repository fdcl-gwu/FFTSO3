#include <omp.h>
#include <iostream>
#include <vector>

#include "fdcl_tictoc.hpp"
#include "fdcl_FFTS2.hpp"


using namespace std;

int main()
{
    int l_max=1000;
    double dt1, dt2;
    fdcl_tictoc tt;
    fdcl_FFTS2_complex FFTS2(l_max);

    omp_set_dynamic(0);     // Explicitly disable dynamic teams

    std::vector<double> dt;
    dt.resize(4);
    for(int r=0; r<=3; r++)
    {
        omp_set_num_threads(std::pow(2,r)); // 
        tt.tic();
        FFTS2.nor_assoc_Legendre_poly(0.1,l_max);
        // FFTS2.spherical_harmonics(0.1,0.2,l_max);
        // FFTS2.compute_weight();
        dt[r]=tt.toc();
    }
    cout << "speed up factor " << dt[0]/dt[1] << " " << dt[0]/dt[2] << endl;

}
