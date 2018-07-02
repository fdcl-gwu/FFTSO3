#include <iostream>
#include <vector>
#include <math.h> // pow
#include <iomanip>      // std::setprecision
#include <Eigen/Dense>
#include <string>

#include "misc_matrix_func.h"

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

        double& operator()(int l, int m, int l1, int m1, int l2, int m2); 
        Eigen::Matrix<double,Dynamic,Dynamic> C;

    private:
        int row(int l, int m, int l1, int m1, int l2, int m2);
        int col(int l, int m, int l1, int m1, int l2, int m2);
};

fdcl_Clebsch_Gordon_matrix::fdcl_Clebsch_Gordon_matrix(int l1, int l2)
{
    init(l1,l2);
}

void fdcl_Clebsch_Gordon_matrix::init(int l1, int l2)
{
    this->l1=l1;
    this->l2=l2;
    C.resize(2*l1+1,2*l2+1);
    C.setZero();
}



int main()
{
    fdcl_Clebsch_Gordon_matrix C(1,2);

    cout << C.l1 << endl;
    return 0;
}

