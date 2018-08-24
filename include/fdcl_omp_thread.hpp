#ifndef _FDCL_OMP_THREAD_HPP
#define _FDCL_OMP_THREAD_HPP

#include <iostream>
#include <omp.h>

namespace fdcl
{
    class omp_thread;
}

class fdcl::omp_thread
{
    public:
        omp_thread(int id, int N_threads);
        ~omp_thread(){};
        int id, N_threads;
        int i_init, i_term, i_init_global, i_term_global;
        void range_open(int i_init_global, int i_term_global);
        void range_closed(int i_init_global, int i_term_global);
};

fdcl::omp_thread::omp_thread(int id, int N_threads)
{
    this->id = id;
    this->N_threads = N_threads;
}

void fdcl::omp_thread::range_open(int i_init_global, int i_term_global)
{
    this->i_init_global = i_init_global;
    this->i_term_global = i_term_global;

    int i_range_global = i_term_global - i_init_global;
    i_init = i_init_global + i_range_global/N_threads * id;
    i_term = i_init + i_range_global/N_threads;
    if (id == N_threads-1)
        i_term = i_term_global;
}

void fdcl::omp_thread::range_closed(int i_init_global, int i_term_global)
{
    this->i_init_global = i_init_global;
    this->i_term_global = i_term_global;

    int i_range_global = i_term_global - i_init_global;
    i_init = i_init_global + i_range_global/N_threads * id;
    i_term = i_init + i_range_global/N_threads-1;
    if (id == N_threads-1)
        i_term = i_term_global;
}
std::ostream& operator<<(std::ostream& os, const fdcl::omp_thread& thr)
{
    os << "thread: " << thr.id << "/" << thr.N_threads-1 << ": [" << thr.i_init << "-" << thr.i_term << "]/[" << thr.i_init_global << "-" << thr.i_term_global << "]";
    return os;
}

#endif
