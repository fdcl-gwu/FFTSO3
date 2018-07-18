#ifndef _FDCL_CLEBSCH_GORDON_MATRIX_HPP
#define _FDCL_CLEBSCH_GORDON_MATRIX_HPP

#include <iostream>
#include <vector>
#include <math.h> // pow
#include <Eigen/Dense>

#include "fdcl_FFTSO3_matrix.hpp"

using namespace std;
using namespace Eigen;

typedef Eigen::Matrix<double, 3, 3> Matrix3;

class fdcl_Clebsch_Gordon_matrix
{
    public:
        fdcl_Clebsch_Gordon_matrix(){};
        fdcl_Clebsch_Gordon_matrix(int l1, int l2);
        void init(int l1, int l2);
        ~fdcl_Clebsch_Gordon_matrix(){};
        int l1, l2;
        Eigen::Matrix<double,Dynamic,Dynamic> C, c;

        double& operator()(int l, int m, int l1, int m1, int l2, int m2); 

        int row(int l, int m, int l1, int m1, int l2, int m2);
        int col(int l, int m, int l1, int m1, int l2, int m2);
        void assert_index(int l, int m, int l1, int m1, int l2, int m2);
        void compute_sub(int l, int m, int l1, int l2);
		void compute(int l1, int l2);
    private:
        fdcl_FFTSO3_matrix_complex matrix2rsph(int );
};

#endif
