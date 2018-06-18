#include "fdcl_FFTSO3_matrix.hpp"

template <class ScalarType>
fdcl_FFTSO3_matrix<ScalarType>::fdcl_FFTSO3_matrix(int l_max)
{
	init(l_max);
}

template <class ScalarType>
void fdcl_FFTSO3_matrix<ScalarType>::init(int l_max)
{
	this->l_max=l_max;
	
	M.resize(l_max+1);
	for(int i=0;i<=l_max;i++)
	{
		M[i].resize(2*i+1,2*i+1);
		M[i].setZero();
	}
}

template <class ScalarType>
void fdcl_FFTSO3_matrix<ScalarType>::assert_index(int l)
{
	assert(l>=0 && l<=l_max);
}

template <class ScalarType>
void fdcl_FFTSO3_matrix<ScalarType>::assert_index(int l, int m, int n)
{
	assert_index(l);
	assert(min(m,n) >= -l && max(m,n) <= l);
}

template <class ScalarType>
Matrix<ScalarType,Dynamic,Dynamic> & fdcl_FFTSO3_matrix<ScalarType>::operator[](int l)
{
	assert_index(l);
	return M[l];
}

template <class ScalarType>
ScalarType& fdcl_FFTSO3_matrix<ScalarType>::operator()(int l, int m, int n)
{
	assert_index(l,m,n);
	return M[l](m+l,n+l);
}

template <class ScalarType>
ostream& operator<< (ostream& os, const fdcl_FFTSO3_matrix<ScalarType>& M)
{
	for(int l=0;l<=M.l_max;l++)
	{
		os << "l=" << l << endl;
		os << M.M[l] << endl << endl;
	}
	return os;
}

/*template<class ScalarType> template<typename U>
fdcl_FFTSO3_matrix<ScalarType> fdcl_FFTSO3_matrix<ScalarType>::operator+(const fdcl_FFTSO3_matrix<U>& M) 
{
	fdcl_FFTSO3_matrix<ScalarType> Z;
	
	Z.init(l_max);
	for(int l=0;l<l_max;l++)
		Z[l]=this->M[l]+M.M[l];
	
	return Z;
}*/

template<class ScalarType>
fdcl_FFTSO3_matrix<complex<double>> fdcl_FFTSO3_matrix<ScalarType>::operator+(fdcl_FFTSO3_matrix<complex<double>> const& M)
{
	fdcl_FFTSO3_matrix<complex<double>> Z(l_max);
	
	for(int l=0;l<=l_max;l++)
		Z[l]=this->M[l]+M.M[l];
			
	return Z;
}

template<class ScalarType>
fdcl_FFTSO3_matrix<ScalarType> fdcl_FFTSO3_matrix<ScalarType>::operator+(fdcl_FFTSO3_matrix<double> const& M)
{
	fdcl_FFTSO3_matrix<ScalarType> Z(l_max);
	
	for(int l=0;l<=l_max;l++)
		Z[l]=this->M[l]+M.M[l];
			
	return Z;
}


template<class ScalarType>
fdcl_FFTSO3_matrix<complex<double>> fdcl_FFTSO3_matrix<ScalarType>::operator*(const complex<double>& c)
{
	fdcl_FFTSO3_matrix<complex<double>> Z;
	
	Z.init(l_max);
	for(int l=0;l<=l_max;l++)
		Z[l]=this->M[l]*c;
			
	return Z;
}

template<class ScalarType>
fdcl_FFTSO3_matrix<ScalarType> fdcl_FFTSO3_matrix<ScalarType>::operator*(const double& c)  	
{
	fdcl_FFTSO3_matrix<ScalarType> Z;
	
	Z.init(l_max);
	for(int l=0;l<=l_max;l++)
		Z[l]=this->M[l]*c;
			
	return Z;
}


template <class ScalarType>
void fdcl_FFTSO3_matrix<ScalarType>::setRandom()
{
	for(int i=0;i<=l_max;i++)
	{
		M[i].setRandom();
	}
}

template <class ScalarType>
void fdcl_FFTSO3_matrix<ScalarType>::setZero()
{
	for(int i=0;i<=l_max;i++)
	{
		M[i].setZero();
	}
}


template class fdcl_FFTSO3_matrix<double>;
template class fdcl_FFTSO3_matrix<complex<double>>;
template ostream& operator<< (ostream& os, const fdcl_FFTSO3_matrix<double>& M);
template ostream& operator<< (ostream& os, const fdcl_FFTSO3_matrix<complex<double>>& M);
//template fdcl_FFTSO3_matrix<double> fdcl_FFTSO3_matrix<double>::operator+ (const fdcl_FFTSO3_matrix<double>& M);
//template fdcl_FFTSO3_matrix<complex<double>> fdcl_FFTSO3_matrix<complex<double>>::operator+ (const fdcl_FFTSO3_matrix<complex<double>>& M);

