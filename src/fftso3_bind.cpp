#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(fftso3, m) {
    m.doc() = "pybind11 example plugin";
}