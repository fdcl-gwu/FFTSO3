cmake_minimum_required(VERSION 3.5)
project( fdcl_FFTSO3 C CXX )

set(CMAKE_CXX_FLAGS "-O3 -Wall -pedantic" )
# set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -pedantic")
# set(CMAKE_CXX_FLAGS_DEBUG "-g")
# Add -Weffc++ for all sorts of good errors to fix
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
#target_compile_definitions(fftso3 PUBLIC EIGEN_FFTW_DEFAULT)
target_link_libraries(fftso3 OpenMP::OpenMP_CXX)
target_link_libraries(fftso3 ${OpenMP_CXX_FLAGS})
#target_link_libraries(fftso3 fftw3 fftw3_threads)

add_executable(test_FFTSO3 src/test_FFTSO3.cpp)
target_link_libraries(test_FFTSO3 fftso3)

add_executable(test_openmp scratch/test_openmp.cpp src/fdcl_tictoc.cpp)   
#target_link_libraries(test_openmp OpenMP::OpenMP_CXX)
#target_link_libraries(test_openmp ${OpenMP_CXX_FLAGS})
target_link_libraries(test_openmp fftso3)

###############################################################################
# GOOGLE TESTING
################################################################################
#
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
    ${CMAKE_BINARY_DIR}/googletest-build
    EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
# if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
# endif()

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(fftso3_unit_test test/test.cpp)
target_link_libraries(fftso3_unit_test gtest_main fftso3)
enable_testing()
add_test(unit_test ../bin/fftso3_unit_test)

