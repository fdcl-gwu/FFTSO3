#include <iostream>
#include <vector>
#include <math.h> // pow, remquo
#include <iomanip>      // std::setprecision
#include <Eigen/Dense>
#include <string> // for class ETOPO5
#include <fstream>

#include "fdcl_FFTSO3.hpp"
#include "fdcl_FFTS2.hpp"
#include "fdcl_tictoc.hpp"
#include "misc_matrix_func.h"

using namespace std;
using namespace Eigen;

typedef Eigen::Matrix<double, 3, 3> Matrix3;

double myrf(double a, double b, double g)
{
    return Euler3232R(a,b,g).trace();
}

double myrfR(Matrix3 R)
{
    return R.trace();
}

complex<double> myf(double a, double b, double g)
{
    return Euler3232R(a,b,g).trace();
}

complex<double> myfR(Matrix3 R)
{
    return R.trace();
}

complex<double> myf_S2(double theta, double phi)
{
    fdcl_FFTS2_complex FFTS2;
    FFTS2.spherical_harmonics(theta,phi,3);

    return cos(theta)*cos(phi)+cos(theta)*sin(phi);
    // return FFTS2.Y(1,1);
    // return 1.0;
}
double myrf_S2(double theta, double phi)
{
    fdcl_FFTS2_complex FFTS2;
    FFTS2.spherical_harmonics(theta,phi,3);

    return cos(theta)*cos(phi)+cos(theta)*sin(phi);
    // return FFTS2.Y(1,1);
    // return 1.0;
}

namespace fdcl
{
    class ETOPO5;
    class spherical_shape_matching;
}

class fdcl::ETOPO5
{
    public:
        ETOPO5(){};
        ~ETOPO5(){};
    
        std::ifstream fd;
        int N_lat, N_lon;
        MatrixXi elev_data;
        VectorXd lat, lon;
        void init(string filename, int N_lat, int N_lon);
        void read(string filename, int N_lat, int N_lon);
        double elev(double lat, double lon);
        double operator()(double lat, double lon); 
    private: 
        void unit_quorem(double num, double den, int& quo, double& rem); 
};

void fdcl::ETOPO5::init(string filename, int N_lat, int N_lon)
{
    this->N_lat=N_lat;
    this->N_lon=N_lon;
    fd.open(filename.c_str());
    elev_data.resize(N_lat,N_lon);
    lat.resize(N_lat);
    lon.resize(N_lon);

    for(int i=0; i<N_lon; i++)
        lon(i)=((double)4320/N_lon)*((double)5/60)*(double)i;
    for(int i=0; i<N_lat; i++)
        lat(i)=((double)2160/N_lat)*((double)5/60)*(double)i;


}

void fdcl::ETOPO5::read(string filename, int N_lat, int N_lon)
{
    init(filename,N_lat,N_lon);
    for(int i_lat=0; i_lat<N_lat; i_lat++)
        for(int i_lon=0; i_lon<N_lon; i_lon++)
            fd >> elev_data(i_lat,i_lon);

}

void fdcl::ETOPO5::unit_quorem(double num, double den, int& quo, double& rem)
{
    double ratio;
    ratio=num/den;
    
    quo=std::floor(ratio);
    rem=ratio-quo;
}

double fdcl::ETOPO5::elev(double lat, double lon)
{
    // bilinear interpolation with scaling
    // https://en.wikipedia.org/wiki/Bilinear_interpolation
    
    assert(lat >= -90. && lat <=90.);
    assert(lon >= 0. && lon <= 360.);

    double x, y, ratio;
    int i0, i1, j0, j1;

    unit_quorem(-lat+90., (((double)2160/N_lat)*((double)5/60)), i0, x);
    unit_quorem(lon, (((double)4320/N_lon)*((double)5/60)), j0, y);

    if (i0==N_lat-1)
        i1=0;
    else
        i1=i0+1;

    if (j0==N_lon-1)
        j1=0;
    else
        j1=j0+1;

    return elev_data(i0,j0)*(1.-x)*(1.-y)+elev_data(i1,j0)*x*(1.-y)+elev_data(i0,j1)*(1.-x)*y+elev_data(i1,j1)*x*y;
}

double fdcl::ETOPO5::operator()(double lat, double lon)
{
    return elev(lat,lon);
}

class fdcl::spherical_shape_matching
{
    public:
        fdcl_FFTS2_matrix_real F, G;
        fdcl_FFTSO3_real RFFTSO3;
        fdcl_FFTSO3_matrix_real U;
        std::vector<fdcl_FFTSO3_matrix_real> u;
        int l_max;
        void init(int l_max);
        double scale_factor=1.e-6, eps=1.e-6;

        spherical_shape_matching(){};
        ~spherical_shape_matching(){};
        spherical_shape_matching(int l_max);

        double J(Eigen::Matrix3d R);
        Eigen::Vector3d gradient(Eigen::Matrix3d R);
        Eigen::Matrix3d opt(Eigen::Matrix3d R0);
        void check_gradient();
};

fdcl::spherical_shape_matching::spherical_shape_matching(int l_max)
{
    init(l_max);
}
void fdcl::spherical_shape_matching::init(int l_max)
{
    this->l_max = l_max;
    F.init(l_max);
    G.init(l_max);
    RFFTSO3.init(l_max);
    U.init(l_max);
    u.resize(4);
    u[1].init(l_max);
    u[2].init(l_max);
    u[3].init(l_max);   
    u = RFFTSO3.deriv_real_harmonics();
 }
double fdcl::spherical_shape_matching::J(Eigen::Matrix3d R)
{
    double y=0.;
    U = RFFTSO3.real_harmonics(R);
    for(int l=0; l<=l_max; l++)
        y+=(G[l]*F[l].transpose()*U[l]).trace();

    return y*scale_factor;
}
Eigen::Vector3d fdcl::spherical_shape_matching::gradient(Eigen::Matrix3d R)
{
    Eigen::Vector3d dJ;
    dJ.setZero();
    U = RFFTSO3.real_harmonics(R);

    for(int l=0; l<=l_max; l++)
        for(int i=1; i<=3; i++)
            dJ(i-1)+=(G[l]*F[l].transpose()*U[l]*u[i][l]).trace();

    return dJ*scale_factor;
}
Eigen::Matrix3d fdcl::spherical_shape_matching::opt(Eigen::Matrix3d R0)
{
    Eigen::Matrix3d R;
    Eigen::Vector3d eta;
    double norm_gradient=1.;
    int i_iter=0;
    R=R0;
    std::vector<double> abg;

    while(norm_gradient > eps)
    {

        eta=gradient(R);
        R=R*expm_SO3(0.0001*eta);
        norm_gradient=eta.norm();
        i_iter+=1;

        abg=R2Euler323(R);

        cout << "iter = " << i_iter << " : J = " << J(R) << " : n_grad = " << norm_gradient << " : Euler = " << abg[0] << " " << abg[1] << " " << abg[2] << endl;
    }

    return R;
}
void fdcl::spherical_shape_matching::check_gradient()
{
    Eigen::Matrix3d R, R_new;
    Eigen::Vector3d eta;
    double a,b,g;

	a=(double)rand()/RAND_MAX*2.*M_PI;
	b=(double)rand()/RAND_MAX*M_PI;
	g=(double)rand()/RAND_MAX*2.*M_PI;

    R=Euler3232R(a,b,g);
    eta.setRandom();
    eta*=1.e-6;

    R_new=R*expm_SO3(eta);

    gradient(R);
    cout << J(R_new)-J(R) << endl;
    cout << gradient(R).transpose()*eta << endl;
}


    
    

    
int main()
{
    int l_max=16;
    fdcl_FFTSO3_complex FFTSO3(l_max);
    fdcl_FFTSO3_real RFFTSO3(l_max);
    fdcl_FFTS2_complex FFTS2(l_max);
    fdcl_FFTS2_real RFFTS2(l_max);
    std::ofstream fd;

    fdcl_tictoc tt;
    double a, b, g;
    a=.12345;
    b=-0.234235;
    g=0.4324235;
    // Matrix3 R;
    // R.setIdentity();
    // Vector3 eta1, eta2;
    // Matrix3 R1, R2;
    // eta1.setRandom();
    // eta2.setRandom();
    // 
    // R1=expm_SO3(eta1);
    // R2=expm_SO3(eta2);

    // FFTSO3.check_verbose=true;
    // FFTSO3.check_all();    
    // RFFTSO3.check_verbose=true;
    // RFFTSO3.check_all();g
    // FFTS2.check_all();
    // RFFTS2.check_all();

    fdcl::ETOPO5 ETOPO5;
    ETOPO5.read("../scratch/ETOPO5.asc",180,360);

    fdcl::spherical_shape_matching SSM(l_max);

    auto func= std::bind(&fdcl::ETOPO5::elev, &ETOPO5, std::placeholders::_1, std::placeholders::_2);
    SSM.F=RFFTS2.forward_transform([=] (double theta, double phi) {return func(-theta*180./M_PI+90.,phi*180./M_PI);} );
    // fd.open("FFTS2.dat");
    // for(int i=0; i<180; i++)
    // {
        // for(int j=0; j<360; j++)
            // fd << RFFTS2.inverse_transform((double)i*M_PI/180., (double)j*M_PI/180.) << " ";
// 
        // fd << endl;
        // cout << i << endl;
// 
    // }
    // fd.close();
// 
    fdcl_FFTSO3_matrix_real U(l_max);
    Eigen::Matrix3d R;

    R=Euler3232R(M_PI/6., -M_PI/3., M_PI/2.);
    U=RFFTSO3.real_harmonics(R);
    for(int l=0; l<=l_max; l++)
        SSM.G[l]=U[l].transpose()*SSM.F[l];

    cout << SSM.J(R) << endl;
    cout << SSM.J(Euler3232R(0.,0.,0.)) << endl;
    SSM.check_gradient();
    cout << SSM.opt(Euler3232R(M_PI/2.,0.,0.)) << endl;
    cout << R << endl;

    // fd.open("FFTS2_rot.dat");
    // for(int i=0; i<180; i++)
    // {
        // for(int j=0; j<360; j++)
            // fd << RFFTS2.inverse_transform((double)i*M_PI/180., (double)j*M_PI/180.) << " ";
// 
        // fd << endl;
        // cout << i << endl;
// 
    // }
    // fd.close();
// 




 

}
