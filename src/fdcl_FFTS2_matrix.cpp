#include "fdcl_FFTS2_matrix.hpp"

template <class ScalarType>
fdcl::FFTS2_matrix<ScalarType>::FFTS2_matrix(int l_max)
{
	init(l_max);
}

template <class ScalarType>
void fdcl::FFTS2_matrix<ScalarType>::init(int l_max)
{
	this->l_max=l_max;
	
	M.resize(l_max+1);
	for(int i=0;i<=l_max;i++)
	{
		M[i].resize(2*i+1);
		M[i].setZero();
	}
}

template <class ScalarType>
void fdcl::FFTS2_matrix<ScalarType>::assert_index(int l)
{
	assert(l>=0 && l<=l_max);
}

template <class ScalarType>
void fdcl::FFTS2_matrix<ScalarType>::assert_index(int l, int m)
{
	assert_index(l);
	assert(m >= -l && m <= l);
}

template <class ScalarType>
Eigen::Matrix<ScalarType,Eigen::Dynamic,1> & fdcl::FFTS2_matrix<ScalarType>::operator[](int l)
{
	assert_index(l);
	return M[l];
}

template <class ScalarType>
ScalarType& fdcl::FFTS2_matrix<ScalarType>::operator()(int l, int m)
{
	assert_index(l,m);
	return M[l](m+l);
}

template <class ScalarType>
ostream& operator<< (ostream& os, const fdcl::FFTS2_matrix<ScalarType>& M)
{
	for(int l=0;l<=M.l_max;l++)
	{
		os << "l=" << l << endl;
		os << M.M[l] << endl << endl;
	}
	return os;
}


template<class ScalarType>
fdcl::FFTS2_matrix<double> fdcl::FFTS2_matrix<ScalarType>::real()  	
{
    fdcl::FFTS2_matrix<double> Z(l_max);

    for(int l=0;l<=l_max;l++)
        Z[l]=this->M[l].real();

    return Z;
}

template<class ScalarType>
fdcl::FFTS2_matrix<complex<double>> fdcl::FFTS2_matrix<ScalarType>::operator+(fdcl::FFTS2_matrix<complex<double>> const& M)
{
	fdcl::FFTS2_matrix<complex<double>> Z(l_max);
	
	for(int l=0;l<=l_max;l++)
		Z[l]=this->M[l]+M.M[l];
			
	return Z;
}

template<class ScalarType>
fdcl::FFTS2_matrix<ScalarType> fdcl::FFTS2_matrix<ScalarType>::operator+(fdcl::FFTS2_matrix<double> const& M)
{
	fdcl::FFTS2_matrix<ScalarType> Z(l_max);
	
	for(int l=0;l<=l_max;l++)
		Z[l]=this->M[l]+M.M[l];
			
	return Z;
}

template<class ScalarType>
fdcl::FFTS2_matrix<complex<double>> fdcl::FFTS2_matrix<ScalarType>::operator-(fdcl::FFTS2_matrix<complex<double>> const& M)
{
	fdcl::FFTS2_matrix<complex<double>> Z(l_max);
	
	for(int l=0;l<=l_max;l++)
		Z[l]=this->M[l]-M.M[l];
			
	return Z;
}

template<class ScalarType>
fdcl::FFTS2_matrix<ScalarType> fdcl::FFTS2_matrix<ScalarType>::operator-(fdcl::FFTS2_matrix<double> const& M)
{
	fdcl::FFTS2_matrix<ScalarType> Z(l_max);
	
	for(int l=0;l<=l_max;l++)
		Z[l]=this->M[l]-M.M[l];
			
	return Z;
}

template<class ScalarType>
fdcl::FFTS2_matrix<complex<double>> fdcl::FFTS2_matrix<ScalarType>::operator*(const complex<double>& c)
{
	fdcl::FFTS2_matrix<complex<double>> Z;
	
	Z.init(l_max);
	for(int l=0;l<=l_max;l++)
		Z[l]=this->M[l]*c;
			
	return Z;
}

template<class ScalarType>
fdcl::FFTS2_matrix<ScalarType> fdcl::FFTS2_matrix<ScalarType>::operator*(const double& c)  	
{
	fdcl::FFTS2_matrix<ScalarType> Z;
	
	Z.init(l_max);
	for(int l=0;l<=l_max;l++)
		Z[l]=this->M[l]*c;
			
	return Z;
}


template <class ScalarType>
void fdcl::FFTS2_matrix<ScalarType>::setRandom()
{
	for(int i=0;i<=l_max;i++)
	{
		M[i].setRandom();
	}
}

template <class ScalarType>
void fdcl::FFTS2_matrix<ScalarType>::setZero()
{
	for(int i=0;i<=l_max;i++)
	{
		M[i].setZero();
	}
}


template <class ScalarType>
double fdcl::FFTS2_matrix<ScalarType>::norm()
{
    double y=0.;

    for(int i=0;i<=l_max;i++)
        y+=M[i].norm();

    return y;
}

namespace fdcl
{
    template class FFTS2_matrix<double>;
    template class FFTS2_matrix<complex<double>>;
}
template ostream& operator<< (ostream& os, const fdcl::FFTS2_matrix<double>& M);
template ostream& operator<< (ostream& os, const fdcl::FFTS2_matrix<complex<double>>& M);

