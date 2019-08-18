namespace fdcl
{
    class spherical_shape_matching;
}

/** Class for Spherical Shape Matching with Harmonic Analysis 
 */
class fdcl::spherical_shape_matching
{
    public:
        fdcl::FFTSO3_real RFFTSO3;

        fdcl::FFTS2_matrix_real F, G; /**< Spherical harmonic for the original and the rotated shape/image */
        fdcl::FFTSO3_matrix_real U; /**< Real harmonics to represent rotation */
        std::vector<fdcl::FFTSO3_matrix_real> u; /**< Derivatives of real harmonics */
        int l_max;
        double scale_factor=1.e-6; /**< scale factor for cost and gradient */
        double eps=1.e-6; /**< stopping threshold */
        double step_size=5.e-3; /**< step size for update */

        spherical_shape_matching(){};
        ~spherical_shape_matching(){f_out.close();};
        spherical_shape_matching(int l_max);

        void init(int l_max); /**< Initialize variables */
        double J(Eigen::Matrix3d R); /**< Cost function for a given rotation */
        Eigen::Vector3d gradient(Eigen::Matrix3d R); /**< Gradient of cost */
        Eigen::Matrix3d opt(Eigen::Matrix3d R0); /**< Steepest descent optimization */
        void check_gradient(); /**< Verify gradient computation */

        std::ofstream f_out;
};

fdcl::spherical_shape_matching::spherical_shape_matching(int l_max)
{
    init(l_max);
    f_out.open("example3_opt.dat");
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
        y+=G[l].transpose()*U[l]*F[l];

    return y*scale_factor;
}

Eigen::Vector3d fdcl::spherical_shape_matching::gradient(Eigen::Matrix3d R)
{
    Eigen::Vector3d dJ;
    dJ.setZero();
    U = RFFTSO3.real_harmonics(R);

    for(int l=0; l<=l_max; l++)
        for(int i=1; i<=3; i++)
            dJ(i-1)+=G[l].transpose()*U[l]*u[i][l]*F[l];

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
        R=R*expm_SO3(step_size*eta);
        norm_gradient=eta.norm();
        i_iter+=1;

        abg=R2Euler323(R);

        cout << "iter = " << i_iter << " : J = " << J(R) << " : n_grad = " << norm_gradient << " : Euler = " << abg[0] << " " << abg[1] << " " << abg[2] << endl;
        f_out << i_iter << ", " << J(R) << ", " << norm_gradient << ", " << abg[0] << ", " << abg[1] << ", " << abg[2] << endl;
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


