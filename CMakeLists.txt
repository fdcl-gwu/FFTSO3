cmake_minimum_required(VERSION 3.5)
project( fdcl_FFTSO3 C CXX )

set(CMAKE_CXX_FLAGS "-O3 -Wall" )
# set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -pedantic")
# set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3")

################################################################################
# GLOBAL SETUP
################################################################################
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

################################################################################
# OpenMP
################################################################################
include(PatchOpenMPApple.cmake)
find_package(OpenMP REQUIRED)
if(NOT TARGET OpenMP::OpenMP_CXX)
   add_library(OpenMP_TARGET INTERFACE)
   add_library(OpenMP::OpenMP_CXX ALIAS OpenMP_TARGET)
   target_compile_options(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
   find_package(Threads REQUIRED)
   target_link_libraries(OpenMP_TARGET INTERFACE Threads::Threads)
   target_link_libraries(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
endif()

################################################################################
# EIGEN
################################################################################
# Eigen setup
# find_package(Eigen3 3.3 REQUIRED)
# if (NOT Eigen3_FOUND)
#      message(STATUS "This project requires Eigen3 and will not be compiled. See http://eigen.tuxfamily.org")
#     return()
# endif()
# include(${EIGEN3_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/extern/eigen3)

################################################################################
# BUILD TARGETS
################################################################################
include_directories(${PROJECT_SOURCE_DIR}/include)

## BUILD FFTSO3 LIBRARY
set(FFTSO3_src 
    src/fdcl_FFTSO3.cpp
    src/fdcl_FFTS2.cpp
    src/fdcl_Clebsch_Gordon.cpp
    src/fdcl_tictoc.cpp
    src/fdcl_FFTSO3_matrix.cpp
    src/fdcl_FFTS2_matrix.cpp
    src/misc_matrix_func.cpp
    src/fdcl_omp_thread.cpp)

add_library(fftso3 SHARED ${FFTSO3_src})
# uncomment the following line to use fftw instead of kissfft, see http://eigen.tuxfamily.org/index.php?title=EigenFFT
# target_compile_definitions(fftso3 PUBLIC EIGEN_FFTW_DEFAULT)
target_link_libraries(fftso3 OpenMP::OpenMP_CXX)
target_link_libraries(fftso3 ${OpenMP_CXX_FLAGS})

## BUILD EXAMPLES
add_executable(example0 example/example0.cpp)
target_link_libraries(example0 fftso3)

add_executable(example1 example/example1.cpp)
target_link_libraries(example1 fftso3)

add_executable(example2 example/example2.cpp)
target_link_libraries(example2 fftso3)

add_executable(example3 example/example3.cpp)
target_link_libraries(example3 fftso3)

add_executable(example4 example/example4.cpp)
target_link_libraries(example4 fftso3)


###############################################################################
# GOOGLE TESTING
################################################################################
#
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(fftso3_unit_test test/test.cpp)
target_link_libraries(fftso3_unit_test gtest_main fftso3)
enable_testing()
add_test(unit_test ../bin/fftso3_unit_test)

