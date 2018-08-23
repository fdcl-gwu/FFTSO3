#include "fdcl_tictoc.hpp"


void fdcl_tictoc::tic()
{
    t0 = std::chrono::steady_clock::now();  
}
double fdcl_tictoc::toc()
{
    double dt;
    t1 = std::chrono::steady_clock::now();  
    dt = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()/1e6;
    cout << "fdcl_tictoc = " << dt << " sec" << endl;
    return dt;
}
double fdcl_tictoc::toc(string message)
{
    cout << message << ": ";
    return toc();
}
