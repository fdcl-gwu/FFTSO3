#include "fdcl_tictoc.hpp"


void fdcl::tictoc::tic()
{
    t0 = std::chrono::steady_clock::now();  
}
double fdcl::tictoc::toc()
{
    double dt;
    t1 = std::chrono::steady_clock::now();  
    dt = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()/1e6;
    std::cout << "fdcl::tictoc = " << dt << " sec" << std::endl;
    return dt;
}
double fdcl::tictoc::toc(std::string message)
{
    std::cout << message << ": ";
    return toc();
}
