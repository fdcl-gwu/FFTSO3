#ifndef _FDCL_TICTOC_HPP
#define _FDCL_TICTOC_HPP

#include <iostream>
#include <string>
#include <chrono>

namespace fdcl
{
    class tictoc;
}

/** \brief Class for stopwatch */
class fdcl::tictoc
{
public:
    std::chrono::steady_clock::time_point t0, t1;
    tictoc(){};
    void tic();
    double toc();
    double toc(std::string message);
    bool quiet=true;
};

#endif
