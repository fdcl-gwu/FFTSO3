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
# set(ENV{OMP_NUM_THREADS} 4)
# set(ENV{OMP_NESTED} TRUE)
# find_package(OpenMP)
# if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

#find_package(Threads REQUIRED)

# Eigen setup
find_package(Eigen3 3.3.4 REQUIRED)
if (NOT Eigen3_FOUND)
    message(STATUS "This project requires Eigen3 and will not be compiled. See http://eigen.tuxfamily.org")
    return()
endif()
include(${EIGEN3_USE_FILE})

################################################################################
# GOOGLE TESTING
################################################################################
#add_subdirectory(${PROJECT_SOURCE_DIR}/extern/googletest/googletest)
#ENABLE_TESTING()
#include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

################################################################################
# BUILD TARGETS
################################################################################
include_directories(${PROJECT_SOURCE_DIR}/include)

set(FFTSO3_src 
    src/test_FFTSO3.cpp
    src/fdcl_FFTSO3.cpp
    src/fdcl_FFTS2.cpp
    src/fdcl_Clebsch_Gordon.cpp
    src/fdcl_tictoc.cpp
    src/fdcl_FFTSO3_matrix.cpp
    src/fdcl_FFTS2_matrix.cpp
    src/misc_matrix_func.cpp)
add_executable(test_FFTSO3 ${FFTSO3_src})

################################################################################
# Testing example
################################################################################
#ENABLE_TESTING()

# add_executable(test_libigl tests/cpp/test_libigl.cpp tests/cpp/test_all.cpp)
# target_link_libraries(test_libigl gtest igl::core igl::cgal)
# add_test(NAME test_libigl COMMAND test_libigl
# WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

