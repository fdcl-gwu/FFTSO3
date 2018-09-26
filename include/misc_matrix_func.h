/*! \file */ 
#ifndef _MISC_MATRIX_FUNC
#define _MISC_MATRIX_FUNC

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace fdcl 
{
    /** Hat map \f$ \hat{\cdot } : \Re^3 \rightarrow \mathfrak{so}(3) \f$
     *
     *  Defined such that \f$ \hat x = x\times y \f$ for all \f$ x,y\in\Re^3 \f$. 
     */
    Eigen::Matrix3d hat(const Eigen::Vector3d v); 

    /** Vee map \f$ \wedge : \mathfrak{so}(3)\rightarrow \Re^3 \f$
     *
     *  The inverse of hat()
     */
    Eigen::Vector3d vee(const Eigen::Matrix3d V);

    /** \f$ \frac{\sin x}{x} \f$
     *
     * Compute \f$ \frac{\sin x}{x} \f$ without the numerical singularity at \f$ x = 0\f$.
     */
    double sinx_over_x(const double x);

    /** Exponential map \f$ \exp: \Re^3 \rightarrow \mathrm{SO(3)} \f$
     *
     * Matrix exponential
     */
    Eigen::Matrix3d expm_SO3(const Eigen::Vector3d r);


    /** Matrix logarithm \f$ \log: \mathrm{SO(3)}\rightarrow \Re^3 \f$
     *
     * The inverse of expm_SO3()
     */
    Eigen::Vector3d logm_SO3(const Eigen::Matrix3d R);

    /** Assert that R is a rotation matrix 
     */
    bool assert_SO3(Eigen::Matrix3d R,const char *R_name);

    /** Convert a rotation matrix into 3-2-3 Euler angles
     */
    std::vector<double> R2Euler323(const Eigen::Matrix3d R);

    /** Convert a set of 3-2-3 Euler angles into a rotation matrix 
     */
    Eigen::Matrix3d Euler3232R(double, double, double);
}

#endif
