#include <iostream>
#include "fdcl_FFTSO3.hpp"

using std::cout;
using std::endl;

double func(fdcl::FFTSO3_matrix_real F, int l_max, double alpha, double beta, double gamma)
{
    fdcl::FFTSO3_real RFFTSO3(l_max);
    return RFFTSO3.inverse_transform(F,alpha,beta,gamma);
}
double func_trace(Eigen::Matrix3d R)
{
    return R.trace();
}
    
int main()
{
    fdcl::FFTSO3_real FFTSO3; 
    fdcl::FFTSO3_matrix_real F, G;
    fdcl::tictoc tt;

/*    std::vector<int> L_max = {4, 8, 16};*/

    //for(auto l_max : L_max)
    //{
        //FFTSO3.init(l_max);
        //F.init(l_max);
        //G.init(l_max);

        //F.setRandom();
        //G=FFTSO3.forward_transform([=] (double alpha, double beta, double gamma)
                //{
                    //return func(F,l_max, alpha,beta,gamma);
                //}
                //);

        //cout << "l_max = " << l_max << ", error = " << (F-G).norm() << endl;
    //}


    // Benchmark : forward transform
    std::vector<int> L_max = {4, 8, 16, 32, 64};
    std::vector<int> N_threads = {1, 2, 4};
    tt.quiet=true;
    int N_repeat = 2;
    double t_elapsed =0.;
    cout << "Benchmark : forward transform " << endl;
    for(auto l_max : L_max)
    {
        for (auto n_threads : N_threads)
        {
            omp_set_num_threads(n_threads);
            FFTSO3.init(l_max);

            tt.tic();
            for(int i_repeat=0; i_repeat < N_repeat; i_repeat++)
                FFTSO3.forward_transform(func_trace);

            t_elapsed = tt.toc();
            t_elapsed /= (double)N_repeat;

            cout << "l_max = " << l_max << ", threads = "<< n_threads << ", t_elapsed " << t_elapsed << endl;
        }
        cout << endl;
    }

    // Benchmark : inverse transform
    cout << "Benchmark : inverse transform " << endl;
    for(auto l_max : L_max)
    {
        for (auto n_threads : N_threads)
        {
            omp_set_num_threads(n_threads);
            FFTSO3.init(l_max);
            F.init(l_max);
            F.setRandom();

            tt.tic();
            for(int i_repeat=0; i_repeat < N_repeat; i_repeat++)
                FFTSO3.inverse_transform(F,0.,0.,0.);

            t_elapsed = tt.toc();
            t_elapsed /= (double)N_repeat;

            cout << "l_max = " << l_max << ", threads = "<< n_threads << ", t_elapsed " << t_elapsed << endl;
        }
        cout << endl;
    }


    return 0;
}


