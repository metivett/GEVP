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
#include "plateau.hpp"

using namespace std;
using namespace LQCDA;
namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;
namespace fs = boost::filesystem;

std::ostream& operator<<(std::ostream& os, const set_scale_parameters& p)
{
	os
    << "pion fit range = {" << p.pi_fit_range.tmin << ", " << p.pi_fit_range.tmax << "}" << endl
    << "K fit range = {" << p.K_fit_range.tmin << ", " << p.K_fit_range.tmax << "}" << endl
    << "Omega fit range = {" << p.Omega_fit_range.tmin << ", " << p.Omega_fit_range.tmax << "}" << endl
    << "nboot = " << p.nboot << endl
    ;
    return os;
}

set_scale_result set_scale(const std::string& manfile, const set_scale_parameters& params)
{
	set_scale_result result;

	string fout_name = manfile.substr(0, manfile.find_last_of('.')) + ".scale.out";
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
 	while(getline(man_if, buf))
 	{
 		ensembles.push_back(buf);
 	}

 	// Extract masses from each ensemble
    int n_ensembles = ensembles.size();

    out << endl << endl
    << "********************************************************************" << endl
    << "Extracting masses from " << n_ensembles << " ensembles :" << endl;
    for(int ens=0; ens<n_ensembles; ens++)
    {
        out << ensembles[ens] << endl;
    }
    out << "********************************************************************" << endl;

    Sample<Matrix<double>> aMpi(params.nboot, n_ensembles, 1);
    Sample<Matrix<double>> aMK(params.nboot, n_ensembles, 1);
    Sample<Matrix<double>> aMOmega(params.nboot, n_ensembles, 1);
 	for(int ens=0; ens<n_ensembles; ens++)
 	{
        out << endl
            << "********************************************************************" << endl
            << "Ensemble " << ensembles[ens] << endl;

        fs::path manfilepath(manfile);
        fs::path conf_list_file_path = manfilepath.parent_path()/ensembles[ens]/"man_pureqcd.GG";
        string conf_list_file = conf_list_file_path.native();

        /*** Pion ***/
        // Read pion correlators for all configurations
        string pi_h("PP00GG");
        Sample<Matrix<double>> pi_corr_sample;
        readSamples(&pi_corr_sample, conf_list_file, pi_h);
        // Resample
        Sample<Matrix<double>> rs_pi_corr = resample(pi_corr_sample, params.nboot, state);
        // Get plateau
        aMpi(ens, 0) = getPlateau(rs_pi_corr, MeffType::COSH, params.pi_fit_range);

        out << "aMpi = " << aMpi.mean()(ens, 0) << endl;

        /*** Kaon ***/
        // Read K correlators for all configurations
        string K_h("PP02GG");
        Sample<Matrix<double>> K_corr_sample;
        readSamples(&K_corr_sample, conf_list_file, K_h);
        // Resample
        Sample<Matrix<double>> rs_K_corr = resample(K_corr_sample, params.nboot, state);
        // Get plateau
        aMK(ens, 0) = getPlateau(rs_K_corr, MeffType::COSH, params.K_fit_range);

        out << "aMK = " << aMK.mean()(ens, 0) << endl;

        /*** Omega ***/
        // Read Omega correlators for all configurations
        string Omega_h("DELTA22GG");
        Sample<Matrix<double>> Omega_corr_sample;
        readSamples(&Omega_corr_sample, conf_list_file, Omega_h);
        // Resample
        Sample<Matrix<double>> rs_Omega_corr = resample(Omega_corr_sample, params.nboot, state);
        // Get plateau
        aMOmega(ens, 0) = getPlateau(rs_Omega_corr, MeffType::SINH, params.Omega_fit_range);

        out << "aMOmega = " << aMOmega.mean()(ens, 0) << endl;
 	}
    out << "********************************************************************" << endl;

    Sample<Matrix<double>> aMKchi2(params.nboot, n_ensembles, 1);
    FOR_SAMPLE(aMKchi2, s)
    {
        // MKchi^2 = MK^2 - Mpi^2 / 2
        aMKchi2[s] = aMK[s].array()*aMK[s].array() - aMpi[s].array()*aMpi[s].array() / 2.;
    }

    Sample<Matrix<double>> aMpi2_aMOmega2(params.nboot, n_ensembles, 1);
    Sample<Matrix<double>> aMKchi2_aMOmega2(params.nboot, n_ensembles, 1);
    FOR_SAMPLE(aMpi2_aMOmega2, s)
    {
        aMpi2_aMOmega2[s] = aMpi[s].array()*aMpi[s].array() / (aMOmega[s].array()*aMOmega[s].array());
        aMKchi2_aMOmega2[s] = aMKchi2[s].array() / (aMOmega[s].array()*aMOmega[s].array());
    }

    /*** Fit scale from Omega ***/
    XYDataSample<double>* aMOmega_vs_aMpi_aMKchi = new XYDataSample<double>(n_ensembles, 2, 1, params.nboot);

    FOR_SAMPLE(aMOmega, s)
    {
        aMOmega_vs_aMpi_aMKchi->x({}, 0)[s] << aMpi2_aMOmega2[s];
        aMOmega_vs_aMpi_aMKchi->x({}, 1)[s] << aMKchi2_aMOmega2[s];
        aMOmega_vs_aMpi_aMKchi->y({}, 0)[s] << aMOmega[s];
    }
    aMOmega_vs_aMpi_aMKchi->setCovFromSample();

    Models::aMOmega_vs_aMpi_aMKchi * model = new Models::aMOmega_vs_aMpi_aMKchi(params.pi_mass, params.K_mass, params.Omega_mass);

    Chi2Fit<double, MIN::MIGRAD> Fit;
    Fit.options.verbosity = SILENT;

    double chi2_dof = 0.;
    vector<double> init_par_values = {0., 0., 0.};

    Sample<double> rs_a(params.nboot);

    FOR_SAMPLE(aMOmega, s)
    {
        // Fit aMOmega_vs_aMpi_aMKchi
        Fit.setData(aMOmega_vs_aMpi_aMKchi->getData(s));
        Fit.fitAllPoints(true);

        auto fit = Fit.fit(*model, init_par_values);
        chi2_dof += fit.cost() / fit.nDOF();

        // // Store results
        double A = fabs(fit.parameters()[0]);
        double A_fit_err = fabs(fit.errors()[0]);
        double B = fabs(fit.parameters()[1]);
        double B_fit_err = fabs(fit.errors()[1]);
        double C = fabs(fit.parameters()[2]);
        double C_fit_err = fabs(fit.errors()[2]);
        rs_a[s] = A;
    }
    chi2_dof /= params.nboot;

    cout << endl << "********************************************************************" << endl;
    cout << "scale setting MOmega fit chi2/dof = " << chi2_dof << endl;

    result.a = rs_a.mean();
    result.a_err = sqrt(rs_a.variance());

	return result;
}