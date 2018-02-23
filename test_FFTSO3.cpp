#include <iostream>
#include <vector>
#include <math.h> // pow
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class fdcl_FFTSO3_matrix
{
public:
	int l_max;
	std::vector<MatrixXd> M;
	fdcl_FFTSO3_matrix(){};
	~fdcl_FFTSO3_matrix(){};
	fdcl_FFTSO3_matrix(int l_max); 
	void init(int l_max);
	MatrixXd& operator[](int l);
	double& operator()(int l, int m, int n);
private:
	void assert_index(int l);
	void assert_index(int l, int m, int n);
};

fdcl_FFTSO3_matrix::fdcl_FFTSO3_matrix(int l_max)
{
	init(l_max);
}

void fdcl_FFTSO3_matrix::init(int l_max)
{
	this->l_max=l_max;
	
	M.resize(l_max+1);
	for(int i=0;i<=l_max;i++)
	{
		M[i].resize(2*i+1,2*i+1);
		M[i].setZero();
	}
}

void fdcl_FFTSO3_matrix::assert_index(int l)
{
	assert(l>=0 && l<=l_max);
}

void fdcl_FFTSO3_matrix::assert_index(int l, int m, int n)
{
	assert_index(l);
	assert(min(m,n) >= -l && max(m,n) <= l);
}

MatrixXd & fdcl_FFTSO3_matrix::operator[](int l)
{
	assert_index(l);
	return M[l];
}

double& fdcl_FFTSO3_matrix::operator()(int l, int m, int n)
{
	assert_index(l,m,n);
	return M[l](m+l,n+l);
}

	
int main()
{
	fdcl_FFTSO3_matrix d;
	
	d.init(3);
	
	d(1,-1,-1)=1.;
	cout << d[1] << endl;

}
