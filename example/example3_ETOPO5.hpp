namespace fdcl
{
    class ETOPO5;
}

/** Access Earth Topological Data, namely ETOPO5 provided by NOAA at https://www.ngdc.noaa.gov/mgg/global/etopo5.HTML
*/

class fdcl::ETOPO5
{
    public:
        ETOPO5(){};
        ~ETOPO5(){};
    
        std::ifstream fd;
        int N_lat, N_lon; /**< Number of grid for latitude and longitude */
        Eigen::MatrixXi elev_data; /**< buffer */
        Eigen::VectorXd lat, lon; /**< latitude and longitude at grid */
        void init(std::string filename, int N_lat, int N_lon); /**< initialize variables */
        void read(std::string filename, int N_lat, int N_lon); /** read data stored at file */
        double elev(double lat, double lon); /**< compute elevation at the given latitude and longitude */
        double operator()(double lat, double lon);  /**< operator overloading for elev() */
    private: 
        void unit_quorem(double num, double den, int& quo, double& rem);  /** compute quotient and remainder */
};

void fdcl::ETOPO5::init(std::string filename, int N_lat, int N_lon)
{
    this->N_lat=N_lat;
    this->N_lon=N_lon;
    fd.open(filename.c_str());
    if (fd.is_open())
    {
        cout << "File successfully open: " << filename << endl;
    }
    else
    {
        cout << "Error opening file: " << filename << endl;
        exit (EXIT_FAILURE);
    }       
    elev_data.resize(N_lat,N_lon);
    lat.resize(N_lat);
    lon.resize(N_lon);

    for(int i=0; i<N_lon; i++)
        lon(i)=((double)4320/N_lon)*((double)5/60)*(double)i;
    for(int i=0; i<N_lat; i++)
        lat(i)=((double)2160/N_lat)*((double)5/60)*(double)i;
}

void fdcl::ETOPO5::read(std::string filename, int N_lat, int N_lon)
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

    double x, y;
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


