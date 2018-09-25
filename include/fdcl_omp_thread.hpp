#ifndef _FDCL_OMP_THREAD_HPP
#define _FDCL_OMP_THREAD_HPP

#include <iostream>
#include <cstdlib> // std::div_t
#include <omp.h>

namespace fdcl
{
    class omp_thread;
}

/** \brief Class to distribute tasks for OpenMP
 *
 * Set of functions to distribute repeated tasks to mutiple threads almost evenly */
class fdcl::omp_thread
{
    public:
        omp_thread(int id, int N_threads);
        ~omp_thread(){};
        int id, N_threads;
        int i_init, i_term, i_init_global, i_term_global;
        void range_open(int i_init_global, int i_term_global);
        void range_closed(int i_init_global, int i_term_global);
    private:
        std::div_t dv;
};

#endif
