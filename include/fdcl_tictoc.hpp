#ifndef _FDCL_TICTOC_HPP
#define _FDCL_TICTOC_HPP

#include <iostream>
#include <string>
#include <chrono>

using namespace std;

class fdcl_tictoc
{
public:
    std::chrono::steady_clock::time_point t0, t1;
    fdcl_tictoc(){};
    void tic();
    void toc();
    void toc(string message);
};

#endif
