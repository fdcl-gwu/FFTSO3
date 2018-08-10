cmake_minimum_required(VERSION 3.9.2)
project( fdcl_FFTSO3 C CXX )

set(CMAKE_CXX_FLAGS "-O3")
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
find_package(Eigen3 3.3 REQUIRED)
if (NOT Eigen3_FOUND)
     message(STATUS "This project requires Eigen3 and will not be compiled. See http://eigen.tuxfamily.org")
     return()
 endif()
include(${EIGEN3_USE_FILE})
# include_directories(${EIGEN3_ROOT_DIR})

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
    src/misc_matrix_func.cpp)

add_library(fftso3 SHARED ${FFTSO3_src})

add_executable(test_FFTSO3 src/test_FFTSO3.cpp)
target_link_libraries(test_FFTSO3 fftso3)

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
add_executable(test test/test.cpp)
target_link_libraries(test gtest_main fftso3)
add_test(NAME example_test COMMAND test)

