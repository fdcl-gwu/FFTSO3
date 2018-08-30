#ifndef _MISC_MATRIX_FUNC
#define _MISC_MATRIX_FUNC

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace fdcl 
{
    Eigen::Matrix3d hat(const Eigen::Vector3d v);
    Eigen::Vector3d vee(const Eigen::Matrix3d V);
    double sinx_over_x(const double x);
    Eigen::Matrix3d expm_SO3(const Eigen::Vector3d r);
    Eigen::Vector3d logm_SO3(const Eigen::Matrix3d R);
    bool assert_SO3(Eigen::Matrix3d R,const char *R_name);
    void sat(Eigen::Vector3d&, double, double );
    void sat(Eigen::Vector4d&, double, double );
    std::vector<double> R2Euler323(const Eigen::Matrix3d R);
    Eigen::Matrix3d Euler3232R(double, double, double);
}

#endif
