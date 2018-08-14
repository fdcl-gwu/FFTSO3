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

int main()
{
    int l_max=127;
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

    auto func= std::bind(&fdcl::ETOPO5::elev, &ETOPO5, std::placeholders::_1, std::placeholders::_2);
    RFFTS2.forward_transform([=] (double theta, double phi) {return func(-theta*180./M_PI+90.,phi*180./M_PI);} );
    fd.open("FFTS2.dat");
    for(int i=0; i<180; i++)
    {
        for(int j=0; j<360; j++)
            fd << RFFTS2.inverse_transform((double)i*M_PI/180., (double)j*M_PI/180.) << " ";

        fd << endl;
        cout << i << endl;

    }
    fd.close();

    fdcl_FFTSO3_matrix_real U(l_max);
    U=RFFTSO3.wigner_D_real(0.,-M_PI/6.,0.,l_max);
    for(int l=0; l<=l_max; l++)
        RFFTS2.F[l]=U[l].transpose()*RFFTS2.F[l];

    fd.open("FFTS2_rot.dat");
    for(int i=0; i<180; i++)
    {
        for(int j=0; j<360; j++)
            fd << RFFTS2.inverse_transform((double)i*M_PI/180., (double)j*M_PI/180.) << " ";

        fd << endl;
        cout << i << endl;

    }
    fd.close();





 

}
