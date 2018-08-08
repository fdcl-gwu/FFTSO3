#include "fdcl_tictoc.hpp"


void fdcl_tictoc::tic()
{
    t0 = std::chrono::steady_clock::now();  
}
void fdcl_tictoc::toc()
{
    t1 = std::chrono::steady_clock::now();  
    cout << "fdcl_tictoc = " << std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()/1e6 << " sec" << endl;
}
void fdcl_tictoc::toc(string message)
{
    cout << message << ": ";
    toc();
}
