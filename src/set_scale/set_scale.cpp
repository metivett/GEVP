/*
 * set_scale.cpp
 *
 *  Created on: Jul 21, 2014
 *      Author: Thibaut Metivet
 */

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/filesystem.hpp>

#include "set_scale.hpp"
#include "extract.hpp"

using namespace std;
using namespace LQCDA;
namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;
namespace fs = boost::filesystem;

std::ostream &operator<<(std::ostream &os, const set_scale_parameters &p)
{
    os
            << "pion fit range = {" << p.pi_fit_range.tmin << ", " << p.pi_fit_range.tmax << "}" << endl
            << "K fit range = {" << p.K_fit_range.tmin << ", " << p.K_fit_range.tmax << "}" << endl
            << "Omega fit range = {" << p.Omega_fit_range.tmin << ", " << p.Omega_fit_range.tmax << "}" << endl
            << "nboot = " << p.nboot << endl
            ;
    return os;
}

std::ostream &operator<<(std::ostream &os, const set_multi_scales_parameters &p)
{
    os
            << "nboot = " << p.nboot << endl
            ;
    return os;
}

Sample<double> set_scale(const std::string &manfile, const set_scale_parameters &params)
{
    fs::path manfilepath(manfile);
    string fout_name = manfilepath.filename().native() + ".scale.out";
    ofstream fout(fout_name);
    tee_device<ostream, ofstream> tee_out(cout, fout);
    stream<tee_device<ostream, ofstream>> out(tee_out);

    // Print info
    out << endl << endl
        << "********************************************************************" << endl
        << "Setting scale with:" << endl
        << "Manifest file " << manfile << endl
        << params
        << "********************************************************************" << endl;

    // Generate random number generator state
    RandGen rng;
    RandGen::rg_state state;
    rng.getState(state);

    // Get ensembles from manfile
    vector<string> ensembles;
    ifstream man_if(manfile);
    string buf;
    while (getline(man_if, buf))
    {
        ensembles.push_back(buf);
    }

    // Extract masses from each ensemble
    int n_ensembles = ensembles.size();

    out << endl << endl
        << "********************************************************************" << endl
        << "Extracting masses from " << n_ensembles << " ensembles :" << endl;
    for (int ens = 0; ens < n_ensembles; ens++)
    {
        out << ensembles[ens] << endl;
    }
    out << "********************************************************************" << endl;

    Sample<Matrix<double>> aMpi(params.nboot, n_ensembles, 1);
    Sample<Matrix<double>> aMK(params.nboot, n_ensembles, 1);
    Sample<Matrix<double>> aMOmega(params.nboot, n_ensembles, 1);
    for (int ens = 0; ens < n_ensembles; ens++)
    {
        out << endl
            << "********************************************************************" << endl
            << "Ensemble " << ensembles[ens] << endl;

        fs::path conf_list_file_path = manfilepath.parent_path() / ensembles[ens] / (ensembles[ens] + ".list");
        string conf_list_file = conf_list_file_path.native();

        /*** Pion ***/
        extract_parameters pi_extract_parameters;
        pi_extract_parameters.particle = Hadron::pion;
        pi_extract_parameters.smearing = "GG";
        pi_extract_parameters.nboot = params.nboot;
        pi_extract_parameters.fold_correlator = true;
        aMpi(ens, 0) = extract_mass(conf_list_file, pi_extract_parameters, MeffType::COSH, params.pi_fit_range);

        out << "aMpi = " << aMpi.mean()(ens, 0) << " +- " << sqrt(aMpi.variance()(ens, 0)) << endl;

        /*** Kaon ***/
        extract_parameters K_extract_parameters;
        K_extract_parameters.particle = Hadron::kaon;
        K_extract_parameters.smearing = "GG";
        K_extract_parameters.nboot = params.nboot;
        K_extract_parameters.fold_correlator = true;
        aMK(ens, 0) = extract_mass(conf_list_file, K_extract_parameters, MeffType::COSH, params.K_fit_range);

        out << "aMK = " << aMK.mean()(ens, 0) << " +- " << sqrt(aMK.variance()(ens, 0)) << endl;

        /*** Omega ***/
        extract_parameters Omega_extract_parameters;
        Omega_extract_parameters.particle = Hadron::omega;
        Omega_extract_parameters.smearing = "GG";
        Omega_extract_parameters.nboot = params.nboot;
        Omega_extract_parameters.fold_correlator = true;
        aMOmega(ens, 0) = extract_mass(conf_list_file, Omega_extract_parameters, MeffType::LOG, params.Omega_fit_range);

        out << "aMOmega = " << aMOmega.mean()(ens, 0) << " +- " << sqrt(aMOmega.variance()(ens, 0)) << endl;
    }
    out << "********************************************************************" << endl;

    Sample<Matrix<double>> aMKchi2(params.nboot, n_ensembles, 1);
    FOR_SAMPLE(aMKchi2, s)
    {
        // MKchi^2 = MK^2 - Mpi^2 / 2
        aMKchi2[s] = aMK[s].array() * aMK[s].array() - aMpi[s].array() * aMpi[s].array() / 2.;
    }

    Sample<Matrix<double>> aMpi2_aMOmega2(params.nboot, n_ensembles, 1);
    Sample<Matrix<double>> aMKchi2_aMOmega2(params.nboot, n_ensembles, 1);
    FOR_SAMPLE(aMpi2_aMOmega2, s)
    {
        aMpi2_aMOmega2[s] = aMpi[s].array() * aMpi[s].array() / (aMOmega[s].array() * aMOmega[s].array());
        aMKchi2_aMOmega2[s] = aMKchi2[s].array() / (aMOmega[s].array() * aMOmega[s].array());
    }

    /*** Fit scale from Omega ***/
    XYDataSample<double> *aMOmega_vs_aMpi2_aMKchi2 = new XYDataSample<double>(n_ensembles, 2, 1, params.nboot);

    FOR_SAMPLE(aMOmega, s)
    {
        aMOmega_vs_aMpi2_aMKchi2->x({}, 0)[s] << aMpi2_aMOmega2[s];
        aMOmega_vs_aMpi2_aMKchi2->x({}, 1)[s] << aMKchi2_aMOmega2[s];
        aMOmega_vs_aMpi2_aMKchi2->y({}, 0)[s] << aMOmega[s];
    }
    // aMOmega_vs_aMpi2_aMKchi2->setCovFromSample();
    Matrix<double> yycov(n_ensembles, n_ensembles);
    yycov.setZero();
    auto aMOmega_var = aMOmega.variance();
    FOR_MAT_DIAG(yycov, i)
    {
        yycov(i, i) = aMOmega_var(i, 0);
    }
    aMOmega_vs_aMpi2_aMKchi2->yyCov(0, 0) = yycov;
    cout << aMOmega_vs_aMpi2_aMKchi2->yyCov(0, 0) << endl;

    Models::aMOmega_vs_aMpi2_aMKchi2 *model = new Models::aMOmega_vs_aMpi2_aMKchi2(params.pi_mass, params.K_mass, params.Omega_mass);

    Chi2Fit<double, MIN::MIGRAD> Fit;
    Fit.options.verbosity = SILENT;

    double chi2_dof = 0., chi2_dof_mean = 0.;
    vector<double> init_par_values = {0.5, 1., 5.};

    // Fit mean aMOmega_vs_aMpi2_aMKchi2
    auto aMOmega_vs_aMpi2_aMKchi2_mean = aMOmega_vs_aMpi2_aMKchi2->mean();
    Fit.setData(aMOmega_vs_aMpi2_aMKchi2_mean);
    Fit.fitAllPoints(true);

    auto fit_mean = Fit.fit(*model, init_par_values);
    chi2_dof_mean += fit_mean.cost() / fit_mean.nDOF();

    // Store results
    double A_mean = (fit_mean.parameters()[0]);
    double A_fit_err_mean = fabs(fit_mean.errors()[0]);
    double B_mean = (fit_mean.parameters()[1]);
    double B_fit_err_mean = fabs(fit_mean.errors()[1]);
    double C_mean = (fit_mean.parameters()[2]);
    double C_fit_err_mean = fabs(fit_mean.errors()[2]);

    init_par_values = {A_mean, B_mean, C_mean};

    Sample<double> rs_a(params.nboot);
    Sample<double> rs_B(params.nboot), rs_C(params.nboot);

    FOR_SAMPLE(aMOmega, s)
    {
        // Fit aMOmega_vs_aMpi2_aMKchi2
        Fit.setData(aMOmega_vs_aMpi2_aMKchi2->getData(s));
        Fit.fitAllPoints(true);

        auto fit = Fit.fit(*model, init_par_values);
        chi2_dof += fit.cost() / fit.nDOF();

        // // Store results
        double A = (fit.parameters()[0]);
        double A_fit_err = fabs(fit.errors()[0]);
        double B = (fit.parameters()[1]);
        double B_fit_err = fabs(fit.errors()[1]);
        double C = (fit.parameters()[2]);
        double C_fit_err = fabs(fit.errors()[2]);
        rs_a[s] = A;
        rs_B[s] = B;
        rs_C[s] = C;
    }
    chi2_dof /= params.nboot;

    cout << endl << "********************************************************************" << endl;
    cout << "scale setting MOmega fit chi2/dof = " << chi2_dof << endl;

    for (int n = 0; n < aMOmega_vs_aMpi2_aMKchi2->nPoints(); n++)
    {
        cout << aMpi2_aMOmega2.mean()(n, 0) << " ";
        cout << sqrt(aMpi2_aMOmega2.variance()(n, 0)) << " ";
        cout << aMKchi2_aMOmega2.mean()(n, 0) << " ";
        cout << sqrt(aMKchi2_aMOmega2.variance()(n, 0)) << " ";
        cout << aMOmega.mean()(n, 0) << " ";
        cout << sqrt(aMOmega.variance()(n, 0)) << " ";
        cout << endl;
    }

    cout << endl << "Mean Fit parameters:" << endl
         << "A = " << A_mean << " +- " << A_fit_err_mean << endl
         << "B = " << B_mean << " +- " << B_fit_err_mean << endl
         << "C = " << C_mean << " +- " << C_fit_err_mean << endl
         << endl;

    cout << endl << "RS Fit parameters:" << endl
         << "A = " << rs_a.mean() << " +- " << sqrt(rs_a.variance()) << endl
         << "B = " << rs_B.mean() << " +- " << sqrt(rs_B.variance()) << endl
         << "C = " << rs_C.mean() << " +- " << sqrt(rs_C.variance()) << endl
         << endl;

    delete aMOmega_vs_aMpi2_aMKchi2;
    delete model;

    return rs_a;
}

Sample<Matrix<double>> set_multi_scales(const std::vector<std::string> &manifests, const set_multi_scales_parameters &params)
{
    string fout_name = "set_multi_scales.out";
    ofstream fout(fout_name);
    tee_device<ostream, ofstream> tee_out(cout, fout);
    stream<tee_device<ostream, ofstream>> out(tee_out);

    // Print info
    cout << endl << endl
         << "********************************************************************" << endl
         << "Setting scale with:" << endl
         << "Manifest files " << endl;
    FOR_VEC(manifests, m)
    {
        cout << manifests[m] << " (beta = " << params.betas[m] << ")" << endl;
    }
    cout << params
         << "********************************************************************" << endl;

    // Generate random number generator state
    RandGen rng;
    RandGen::rg_state state;
    rng.getState(state);

    // Process all manifest files (extract pion, K, omega plateaus)
    Sample<Matrix<double>> aMpi(params.nboot);
    Sample<Matrix<double>> aMK(params.nboot);
    Sample<Matrix<double>> aMOmega(params.nboot);
    Sample<Matrix<double>> beta_index(params.nboot);

    std::unordered_map<double, unsigned int> beta_indexes;

    unsigned int pt_index = 0;

    for (int man = 0; man < manifests.size(); man++)
    {
        std::string manfile = manifests[man];
        fs::path manfilepath(manfile);
        double cur_beta = params.betas[man];

        if (beta_indexes.count(cur_beta) == 0)
        {
            auto b_ind = beta_indexes.size();
            beta_indexes[cur_beta] = b_ind;
        }

        out << endl << endl
            << "********************************************************************" << endl
            << "********************************************************************" << endl
            << "Processing manifest file " << manfile << endl;

        // Get ensembles from manfile
        vector<string> ensembles;
        ifstream man_if(manfile);
        string buf;
        while (getline(man_if, buf))
        {
            ensembles.push_back(buf);
        }

        // Extract masses from each ensemble
        int n_ensembles = ensembles.size();

        out << endl << endl
            << "********************************************************************" << endl
            << "Extracting masses from " << n_ensembles << " ensembles :" << endl;
        for (int ens = 0; ens < n_ensembles; ens++)
        {
            out << ensembles[ens] << endl;
        }
        out << "********************************************************************" << endl;

        aMpi.conservativeResizeMatrix(aMpi.rows() + n_ensembles, 1);
        aMK.conservativeResizeMatrix(aMK.rows() + n_ensembles, 1);
        aMOmega.conservativeResizeMatrix(aMOmega.rows() + n_ensembles, 1);
        beta_index.conservativeResizeMatrix(beta_index.rows() + n_ensembles, 1);

        for (int ens = 0; ens < n_ensembles; ens++)
        {
            out << endl
                << "********************************************************************" << endl
                << "Ensemble " << ensembles[ens] << endl;

            fs::path conf_list_file_path = manfilepath.parent_path() / ensembles[ens] / (ensembles[ens] + ".list");
            string conf_list_file = conf_list_file_path.native();

            /*** Pion ***/
            extract_parameters pi_extract_parameters;
            pi_extract_parameters.particle = Hadron::pion;
            pi_extract_parameters.smearing = "GG";
            pi_extract_parameters.nboot = params.nboot;
            pi_extract_parameters.beta = cur_beta;
            pi_extract_parameters.fold_correlator = true;
            aMpi(pt_index, 0) = extract_mass(conf_list_file, pi_extract_parameters, MeffType::COSH, fit_range());

            out << "aMpi = " << aMpi.mean()(pt_index, 0) << " +- " << sqrt(aMpi.variance()(pt_index, 0)) << endl;

            /*** Kaon ***/
            extract_parameters K_extract_parameters;
            K_extract_parameters.particle = Hadron::kaon;
            K_extract_parameters.smearing = "GG";
            K_extract_parameters.nboot = params.nboot;
            K_extract_parameters.beta = cur_beta;
            K_extract_parameters.fold_correlator = true;
            aMK(pt_index, 0) = extract_mass(conf_list_file, K_extract_parameters, MeffType::COSH, fit_range());

            out << "aMK = " << aMK.mean()(pt_index, 0) << " +- " << sqrt(aMK.variance()(pt_index, 0)) << endl;

            /*** Omega ***/
            extract_parameters Omega_extract_parameters;
            Omega_extract_parameters.particle = Hadron::omega;
            Omega_extract_parameters.smearing = "GG";
            Omega_extract_parameters.nboot = params.nboot;
            Omega_extract_parameters.beta = cur_beta;
            Omega_extract_parameters.fold_correlator = false;
            aMOmega(pt_index, 0) = extract_mass(conf_list_file, Omega_extract_parameters, MeffType::LOG, fit_range());

            out << "aMOmega = " << aMOmega.mean()(pt_index, 0) << " +- " << sqrt(aMOmega.variance()(pt_index, 0)) << endl;

            /*** beta of the considered point ***/
            FOR_SAMPLE(beta_index, s)
            {
                beta_index[s](pt_index, 0) = beta_indexes[cur_beta];
            }

            pt_index++;

        } // for (int ens = 0; ens < n_ensembles; ens++)

        out << "********************************************************************" << endl;
    } // for (int man = 0; man < manifests.size(); man++)

    unsigned int npoints = pt_index;

    out << endl << endl
        << "********************************************************************" << endl
        << "********************************************************************" << endl;

    Sample<Matrix<double>> aMKchi2(params.nboot, npoints, 1);
    FOR_SAMPLE(aMKchi2, s)
    {
        // MKchi^2 = MK^2 - Mpi^2 / 2
        aMKchi2[s] = aMK[s].array() * aMK[s].array() - aMpi[s].array() * aMpi[s].array() / 2.;
    }

    Sample<Matrix<double>> aMpi2_aMOmega2(params.nboot, npoints, 1);
    Sample<Matrix<double>> aMKchi2_aMOmega2(params.nboot, npoints, 1);
    FOR_SAMPLE(aMpi2_aMOmega2, s)
    {
        aMpi2_aMOmega2[s] = aMpi[s].array() * aMpi[s].array() / (aMOmega[s].array() * aMOmega[s].array());
        aMKchi2_aMOmega2[s] = aMKchi2[s].array() / (aMOmega[s].array() * aMOmega[s].array());
    }

    /*** Fit scale from Omega ***/
    XYDataSample<double> *aMOmega_vs_aMpi2_aMKchi2 = new XYDataSample<double>(npoints, 3, 1, params.nboot);

    FOR_SAMPLE(aMOmega, s)
    {
        aMOmega_vs_aMpi2_aMKchi2->x({}, 0)[s] << aMpi2_aMOmega2[s];
        aMOmega_vs_aMpi2_aMKchi2->x({}, 1)[s] << aMKchi2_aMOmega2[s];
        aMOmega_vs_aMpi2_aMKchi2->x({}, 2)[s] << beta_index[s];
        aMOmega_vs_aMpi2_aMKchi2->y({}, 0)[s] << aMOmega[s];
    }
    // aMOmega_vs_aMpi2_aMKchi2->setCovFromSample();
    Matrix<double> yycov(npoints, npoints);
    yycov.setZero();
    auto aMOmega_var = aMOmega.variance();
    FOR_MAT_DIAG(yycov, i)
    {
        yycov(i, i) = aMOmega_var(i, 0);
    }
    aMOmega_vs_aMpi2_aMKchi2->yyCov(0, 0) = yycov;
    // cout << aMOmega_vs_aMpi2_aMKchi2->yyCov(0, 0) << endl;

    unsigned int nbetas = beta_indexes.size();

    Models::aMOmega_vs_aMpi2_aMKchi2_multi *model = new Models::aMOmega_vs_aMpi2_aMKchi2_multi(params.pi_mass, params.K_mass, params.Omega_mass, nbetas);

    Chi2Fit<double, MIN::MIGRAD> Fit;
    Fit.options.verbosity = SILENT;

    double chi2_dof = 0., chi2_dof_mean = 0.;
    vector<double> init_par_values(nbetas + 2);
    for (int nb = 0; nb < nbetas; nb++)
    {
        init_par_values[nb] = 0.5;
    }
    init_par_values[nbetas] = 1.;
    init_par_values[nbetas + 1] = 6.;


    // Fit mean aMOmega_vs_aMpi2_aMKchi2
    auto aMOmega_vs_aMpi2_aMKchi2_mean = aMOmega_vs_aMpi2_aMKchi2->mean();
    Fit.setData(aMOmega_vs_aMpi2_aMKchi2_mean);
    Fit.fitAllPoints(true);

    auto fit_mean = Fit.fit(*model, init_par_values);
    chi2_dof_mean += fit_mean.cost() / fit_mean.nDOF();

    // Store results
    std::vector<double> A_mean(nbetas), A_fit_err_mean(nbetas);
    for (int nb = 0; nb < nbetas; nb++)
    {
        A_mean[nb] = (fit_mean.parameters()[nb]);
        A_fit_err_mean[nb] = fabs(fit_mean.errors()[nb]);
    }
    double B_mean = (fit_mean.parameters()[nbetas]);
    double B_fit_err_mean = fabs(fit_mean.errors()[nbetas]);
    double C_mean = (fit_mean.parameters()[nbetas + 1]);
    double C_fit_err_mean = fabs(fit_mean.errors()[nbetas + 1]);

    std::copy(A_mean.begin(), A_mean.end(), init_par_values.begin());
    init_par_values[nbetas] = B_mean;
    init_par_values[nbetas + 1] = C_mean;

    Sample<Matrix<double>> rs_a(params.nboot, nbetas, 1);
    Sample<double> rs_B(params.nboot), rs_C(params.nboot);

    FOR_SAMPLE(aMOmega, s)
    {
        // Fit aMOmega_vs_aMpi2_aMKchi2
        Fit.setData(aMOmega_vs_aMpi2_aMKchi2->getData(s));
        Fit.fitAllPoints(true);

        auto fit = Fit.fit(*model, init_par_values);
        chi2_dof += fit.cost() / fit.nDOF();

        // Store results
        std::vector<double> A(nbetas);
        for (int nb = 0; nb < nbetas; nb++)
        {
            A[nb] = (fit.parameters()[nb]);
            double A_fit_err = fabs(fit.errors()[nb]);
        }
        double B = (fit.parameters()[nbetas]);
        double B_fit_err = fabs(fit.errors()[nbetas]);
        double C = (fit.parameters()[nbetas + 1]);
        double C_fit_err = fabs(fit.errors()[nbetas + 1]);
        for (int nb = 0; nb < nbetas; nb++)
            rs_a[s](nb, 0) = A[nb];
        rs_B[s] = B;
        rs_C[s] = C;
    }
    chi2_dof /= params.nboot;

    cout << endl << "********************************************************************" << endl;
    cout << "scale setting MOmega fit chi2/dof = " << chi2_dof << endl;

    for (int n = 0; n < aMOmega_vs_aMpi2_aMKchi2->nPoints(); n++)
    {
        cout << aMpi2_aMOmega2.mean()(n, 0) << " ";
        cout << sqrt(aMpi2_aMOmega2.variance()(n, 0)) << " ";
        cout << aMKchi2_aMOmega2.mean()(n, 0) << " ";
        cout << sqrt(aMKchi2_aMOmega2.variance()(n, 0)) << " ";
        cout << aMOmega.mean()(n, 0) << " ";
        cout << sqrt(aMOmega.variance()(n, 0)) << " ";
        cout << beta_index.mean()(n, 0) << " ";
        cout << endl;
    }

    cout << endl << "Mean Fit parameters:" << endl;
    for (int nb = 0; nb < nbetas; nb++)
    {
        cout << "A" << nb << " = " << A_mean[nb] << " +- " << A_fit_err_mean[nb] << endl;
    }
    cout << "B = " << B_mean << " +- " << B_fit_err_mean << endl
         << "C = " << C_mean << " +- " << C_fit_err_mean << endl
         << endl;

    cout << endl << "RS Fit parameters:" << endl;
    for (int nb = 0; nb < nbetas; nb++)
    {
        cout << "A" << nb << " = " << rs_a.mean()(nb, 0) << " +- " << sqrt(rs_a.variance()(nb, 0)) << endl;
    }
    cout << "B = " << rs_B.mean() << " +- " << sqrt(rs_B.variance()) << endl
         << "C = " << rs_C.mean() << " +- " << sqrt(rs_C.variance()) << endl
         << endl;

    delete aMOmega_vs_aMpi2_aMKchi2;
    delete model;

    return rs_a;
}