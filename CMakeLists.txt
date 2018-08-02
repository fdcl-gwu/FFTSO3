cmake_minimum_required(VERSION 3.9.2)
project( fdcl_FFTSO3 C CXX )

# set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
# Add -Weffc++ for all sorts of good errors to fix
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3")
################################################################################
# GLOBAL SETUP
################################################################################
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Find OpenMP
set(ENV{OMP_NUM_THREADS} 4)
set(ENV{OMP_NESTED} TRUE)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(Threads REQUIRED)

# Eigen setup
find_package(Eigen3 3.3.4 REQUIRED)
if (NOT Eigen3_FOUND)
    message(STATUS "This project requires Eigen3 and will not be compiled.")
    return()
endif()
include(${EIGEN3_USE_FILE})

# just link executables against HDF5 target_link_libraries(target ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

# setup Python
# include(FindPythonAnaconda)
# Pybind11 setup
find_package(PythonInterp REQUIRED)
find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)
if (PYTHONINTERP_FOUND)
    message(STATUS "Python Found: ${PYTHON_VERSION_STRING}")
endif()
include_directories(${PYTHON_INCLUDE_DIRS})
link_libraries(${PYTHON_LIBRARIES})
link_directories(${PYTHON_LIBRARY_DIRS})

add_subdirectory(extern/pybind11)

################################################################################
# GOOGLE TESTING
################################################################################
add_subdirectory(${PROJECT_SOURCE_DIR}/extern/googletest/googletest)
ENABLE_TESTING()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

################################################################################
# BUILD TARGETS
################################################################################
include_directories(${PROJECT_SOURCE_DIR}/include)

set(test_FFTSO3_src 
    src/test_FFTSO3.cpp
    src/fdcl_FFTSO3_real.cpp
    src/fdcl_FFTSO3_complex.cpp
    src/fdcl_Clebsch_Gordon_matrix.cpp
    src/fdcl_tictoc.cpp
    src/fdcl_FFTSO3_matrix.cpp
    src/misc_matrix_func.cpp)
add_executable(test_FFTSO3 ${test_FFTSO3_src})
# mesh remeshing example
# add_executable(remesh src/remesh_main.cpp)
# target_link_libraries(remesh cgal_cpp input_parser_cpp igl::core)
################################################################################
# PYTHON BINDING
################################################################################
pybind11_add_module(fftso3 MODULE src/fftso3_bind.cpp src/misc_matrix_func)
# target_link_libraries(asteroid PRIVATE igl::core cgal_cpp)

################################################################################
# Testing example
################################################################################
ENABLE_TESTING()

# add_executable(test_libigl tests/cpp/test_libigl.cpp tests/cpp/test_all.cpp)
# target_link_libraries(test_libigl gtest igl::core igl::cgal)
# add_test(NAME test_libigl COMMAND test_libigl
# WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

