#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "misc_matrix_func.h"

namespace py = pybind11;


PYBIND11_MODULE(fftso3, m) {
    m.doc() = "pybind11 FFTSO3 plugin";
    m.def("hat", &hat, "hat");
    m.def("vee", &vee);
    m.def("sinx_over_x", &sinx_over_x);
    m.def("expm_SO3", &expm_SO3);
    m.def("logm_SO3", &logm_SO3);
    m.def("assert_SO3", &assert_SO3);
    m.def("R2Euler323", &R2Euler323);
    m.def("Euler3232R", &Euler3232R);
}