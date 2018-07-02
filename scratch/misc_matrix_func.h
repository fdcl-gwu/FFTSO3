#ifndef _MISC_MATRIX_FUNC
#define _MISC_MATRIX_FUNC

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using Eigen::MatrixXd;
using namespace std;

typedef Eigen::Matrix<double, 3, 3> Matrix3;
typedef Eigen::Matrix<double, 3, 1> Vector3;
typedef Eigen::Matrix<double, 4, 1> Vector4;

Matrix3 hat(const Vector3 v);
Vector3 vee(const Matrix3 V);
double sinx_over_x(const double x);
Matrix3 expm_SO3(const Vector3 r);
Vector3 logm_SO3(const Matrix3 R);
bool assert_SO3(Matrix3 R,const char *R_name);
void sat(Vector3&, double, double );
void sat(Vector4&, double, double );
std::vector<double> R2Euler323(const Matrix3 R);
Matrix3 Euler3232R(double, double, double);

#endif
