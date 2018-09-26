#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "misc_matrix_func.h"

Eigen::Matrix3d fdcl::hat(const Eigen::Vector3d v)
{
	Eigen::Matrix3d V;
	V.setZero();
	V(2,1)=v(0);V(1,2)=-V(2,1);
	V(0,2)=v(1);V(2,0)=-V(0,2);
	V(1,0)=v(2);V(0,1)=-V(1,0);

	return V;
}

Eigen::Vector3d fdcl::vee(const Eigen::Matrix3d V)
{
	Eigen::Vector3d v;
	Eigen::Matrix3d E;
	v.setZero();
	E=V+V.transpose();
	if(E.norm() > 1.e-6)
	{
        std::cout << "vee: ERROR: the input matrix is not skew-symmetric" << std::endl;
	}
	else
	{
		v(0)=V(2,1);
		v(1)=V(0,2);
		v(2)=V(1,0);
	}
	return  v;
}

double fdcl::sinx_over_x(const double x)
{
	double y;
	double eps=1.e-6;
	if(abs(x)<eps)
	    y=- pow(x,10)/39916800. + pow(x,8)/362880. - pow(x,6)/5040. + pow(x,4)/120. - pow(x,2)/6. + 1.;
	else
		y=sin(x)/x;

	return y;
}

Eigen::Matrix3d fdcl::expm_SO3(const Eigen::Vector3d r)
{
	Eigen::Matrix3d R;
	double theta,y,y2;

	theta=r.norm();
	y=sinx_over_x(theta);
	y2=sinx_over_x(theta/2.);

	R.setIdentity();
	R+=y*hat(r)+1./2.*pow(y2,2)*hat(r)*hat(r);

	return R;
}

Eigen::Vector3d fdcl::logm_SO3(const Eigen::Matrix3d R)
{
	Eigen::Vector3d r;
	Eigen::Matrix3d I;
	double eps=1.e-6;

	r.setZero();
	I.setIdentity();
	if((I-R*R.transpose()).norm() > eps || abs(R.determinant()-1) > eps)
	{
        std::cout << "logm_SO3: error: R is not a rotation matrix" << std::endl;
	}
	else
	{
		Eigen::EigenSolver<Eigen::MatrixXd> eig(R);
		double min_del_lam_1=1.0, cos_theta, theta;
		int i,i_min=-1;
		Eigen::Vector3d v;
		Eigen::Matrix3d R_new;
		for(i=0;i<3;i++)
		{
			if(eig.eigenvectors().col(i).imag().norm() < eps)
			{
				if(pow(eig.eigenvalues()[i].real(),2)-1. < min_del_lam_1 )
				{
					min_del_lam_1=pow(eig.eigenvalues()[i].real(),2)-1.0;
					i_min=i;
				}
			}

		}
		v=eig.eigenvectors().col(i_min).real();

		cos_theta=(R.trace()-1.0)/2.0;
		if(cos_theta > 1.0)
			cos_theta=1.0;
		else if(cos_theta < -1.0)
			cos_theta=-1.0;
		theta=acos(cos_theta);
		R_new=expm_SO3(theta*v);

		if((R-R_new).norm() > (R-R_new.transpose()).norm())
			v=-v;

		r=theta*v;

		return r;

	}

	return r;
}

bool fdcl::assert_SO3(Eigen::Matrix3d R,const char *R_name)
{
	bool isSO3;
	double errO, errD;
	Eigen::Matrix3d eye3;
	double eps=1e-5;

	eye3.setIdentity();
	errO=(R.transpose()*R-eye3).norm();
	errD=pow(1-R.determinant(),2);



	if (errO > eps || errD > eps)
	{
		isSO3=false;
		printf("ERROR: %s: ||I-R^TR||= %8.6f, det(R)= %8.6f\n",R_name,errO,R.determinant());
        std::cout << R << std::endl << std::endl;
	}
	else
	{
		isSO3=true;
	}
	return isSO3;

}

std::vector<double> fdcl::R2Euler323(const Eigen::Matrix3d R)
{
	double a, b, g;
	std::vector<double> abg;

	b=acos(R(2,2));

	if(b!=0.)
	{
	    g=atan2(R(2,1),-R(2,0));
	    a=atan2(R(1,2),R(0,2));
	}
	else
	{
	    a=acos(R(0,0))/2.;
	    g=a;
	}

	abg.resize(3);
	abg[0]=a;
	abg[1]=b;
	abg[2]=g;
	
	return abg;
}
Eigen::Matrix3d fdcl::Euler3232R(double a, double b, double g)
{
	Eigen::Matrix3d R;
	
	R << cos(a)*cos(b)*cos(g) - sin(a)*sin(g), - cos(g)*sin(a) - cos(a)*cos(b)*sin(g), cos(a)*sin(b),
		cos(a)*sin(g) + cos(b)*cos(g)*sin(a),   cos(a)*cos(g) - cos(b)*sin(a)*sin(g), sin(a)*sin(b),
		-cos(g)*sin(b),                          sin(b)*sin(g),        cos(b);

	return R;
}
