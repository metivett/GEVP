/*
 * analyze.cpp
 *
 *  Created on: Apr 11, 2014
 *      Author: Thibaut Metivet
 */

 #include <iostream>
 #include <complex>
 #include <boost/iostreams/tee.hpp>
 #include <boost/iostreams/stream.hpp>

 #include "LQCDA.hpp"
 #include "analyze.hpp"
 #include "utils.hpp"
 #include "plateau.hpp"
 #include "gevp.hpp"

 using namespace std;
 using namespace LQCDA;
 namespace bio = boost::iostreams;
 using bio::tee_device;
 using bio::stream;

 std::ostream& operator<<(std::ostream& os, const analysis_parameters& p)
 {
 	std::cout
 	<< "L = " << p.L << endl
 	<< "T = " << p.T << endl
 	<< "beta = " << p.beta << endl
 	<< "frame = " << p.frame << endl
 	<< "t0 = " << p.t0 << endl
 	<< "nboot = " << p.nboot << endl
    << "local fit range = {" << p.local_range.tmin << ", " << p.local_range.tmax << "}" << endl
    << "gevp fit range = {" << p.gevp_range.tmin << ", " << p.gevp_range.tmax << "}" << endl;
 }

 int analyze(const std::string& manfile, const analysis_parameters& params)
 {
    string fout_name = manfile.substr(0, manfile.find_last_of('.')) + ".out";
    ofstream fout(fout_name);
    tee_device<ostream, ofstream> tee_out(cout, fout);
    stream<tee_device<ostream, ofstream>> out(tee_out);

 	// Print info
 	out << endl << endl
 	<< "********************************************************************" << endl
 	<< "Performing GEVP analysis with :" << endl
 	<< "Manifest file " << manfile << endl
 	<< params
 	<< "********************************************************************" << endl;

 	// Generate random number generator state
 	RandGen rng;
 	RandGen::rg_state state;
 	rng.getState(state);

    // Time range
    int Nt = params.T / 2;


 	/*** Local pion operator analysis (COM) ***/

 	out << endl
 	<< "********************************************************************" << endl
 	<< "Performing Local pion operator analysis :" << endl;

    // Read PionPion correlators for all configurations
 	std::string pi_h("PP00GG");

    out << "reading samples... ";
 	Sample<Matrix<double>> pi_corr_sample;
 	readSamples(&pi_corr_sample, manfile, pi_h);
    out << "done\n";

 	Matrix<double> pi_t = pi_corr_sample[0].col(0);

    // Resample PP correlators
    out << "resampling... ";
 	Sample<Matrix<double>> rs_pi_corr = resample(pi_corr_sample, params.nboot, state);
    out << "done\n";

    Vector<double> pi_corr = rs_pi_corr.mean().col(0);
    Vector<double> pi_corr_err = rs_pi_corr.variance().col(0);

    // Print correlator
    ofstream of_pi_corr("pi_corr.dat");
    for(int i = 0; i < pi_corr.size(); ++i) {
        of_pi_corr << i << ' ' << pi_corr[i] << ' ' << sqrt(pi_corr_err[i]) << endl;
    }

    // Effective mass
    out << "computing effective mass... ";
 	Sample<Matrix<double>> rs_pi_meff(params.nboot);
 	rs_pi_meff.resizeMatrix(Nt, 1);
 	for(int s = 0; s < params.nboot; ++s) {
 		for(int i = 0; i < Nt; ++i) {
 			rs_pi_meff[s](i) = (std::log(rs_pi_corr[s](i, 0) / rs_pi_corr[s](i+1, 0)));
 		}
 	}

 	Vector<double> pi_meff = rs_pi_meff.mean();
 	Vector<double> pi_meff_err = rs_pi_meff.variance();
    out << "done\n";

    // Print effective mass
 	ofstream of_pi_meff("pi_meff.dat");
 	for(int i = 0; i < pi_meff.size(); ++i) {
 		of_pi_meff << i << ' ' << pi_meff[i] << ' ' << sqrt(pi_meff_err[i]) << endl;
 	}

    // Compute resampled plateaus
    out << "fitting plateau... \n";
 	Sample<double> rs_pi_plat = fitPlateau(rs_pi_meff, params.local_range);
    out << "done\n";

 	double pi_plat = rs_pi_plat.mean();
 	double pi_plat_err = rs_pi_plat.variance();

    // Print results
 	out << "Pion analysis : " << endl
 	<< "E_COM = " << pi_plat << " (" << GeV(pi_plat, params.beta) << " GeV)"
 	<< " +- " << sqrt(pi_plat_err) << " (" << GeV(sqrt(pi_plat_err), params.beta) << " GeV)"
 	<< endl 
    << "********************************************************************" << endl << endl;

 	out << endl;



    /*** Local rho operator analysis (COM) ***/

    out << endl
    << "********************************************************************" << endl
    << "Performing Local rho operator analysis :" << endl;

    // Read V1V1 correlators for all configurations
    std::string rho_h("V1V100GG");

    out << "reading samples... ";
    Sample<Matrix<double>> rho_corr_sample;
    readSamples(&rho_corr_sample, manfile, rho_h);
    out << "done\n";

    Matrix<double> rho_t = rho_corr_sample[0].col(0);

    // Resample V1V1 correlators
    out << "resampling... ";
    Sample<Matrix<double>> rs_rho_corr = resample(rho_corr_sample, params.nboot, state);
    out << "done\n";

    Vector<double> rho_corr = rs_rho_corr.mean().col(0);
    Vector<double> rho_corr_err = rs_rho_corr.variance().col(0);

    // Print correlator
    ofstream of_rho_corr("rho_corr.dat");
    for(int i = 0; i < rho_corr.size(); ++i) {
        of_rho_corr << i << ' ' << rho_corr[i] << ' ' << sqrt(rho_corr_err[i]) << endl;
    }

    // Effective mass
    out << "computing effective mass... ";
    Sample<Matrix<double>> rs_rho_meff(params.nboot);
    rs_rho_meff.resizeMatrix(Nt, 1);
    for(int s = 0; s < params.nboot; ++s) {
        for(int i = 0; i < Nt; ++i) {
            rs_rho_meff[s](i) = (std::log(rs_rho_corr[s](i, 0) / rs_rho_corr[s](i+1, 0)));
        }
    }

    Vector<double> rho_meff = rs_rho_meff.mean();
    Vector<double> rho_meff_err = rs_rho_meff.variance();
    out << "done\n";

    // Print effective mass
    ofstream of_rho_meff("rho_meff.dat");
    for(int i = 0; i < rho_meff.size(); ++i) {
        of_rho_meff << i << ' ' << rho_meff[i] << ' ' << sqrt(rho_meff_err[i]) << endl;
    }

    // Compute resampled plateaus
    out << "fitting plateau... \n";
    Sample<double> rs_rho_plat = fitPlateau(rs_rho_meff, params.local_range);
    out << "done\n";

    double rho_plat = rs_rho_plat.mean();
    double rho_plat_err = rs_rho_plat.variance();

    // Print results
    out << "Rho analysis : " << endl
    << "E_COM = " << rho_plat << " (" << GeV(rho_plat, params.beta) << " GeV)"
    << " +- " << sqrt(rho_plat_err) << " (" << GeV(sqrt(rho_plat_err), params.beta) << " GeV)"
    << endl
    << "********************************************************************" << endl << endl;

    out << endl;


    /*** GEVP analysis (COM) ***/

    out << endl
    << "********************************************************************" << endl
    << "Performing rho-pipi GEVP analysis :" << endl;

    int t0 = params.t0;

    // Read ViVi correlators for all configurations
    std::string gevp_h("ViVi");

    out << "reading samples... ";
    Sample<Matrix<double>> gevp_corr_sample;
    readSamples(&gevp_corr_sample, manfile, gevp_h);
    out << "done\n";

    Matrix<double> gevp_t = gevp_corr_sample[0].col(0);

    int ngev = sqrt((gevp_corr_sample[0].cols() - 1) / 2);

    // Resample ViVi correlators
    out << "resampling... ";
    Sample<Matrix<double>> rs_gevp_corr = resample(gevp_corr_sample, params.nboot, state);
    out << "done\n";

    // Replace (rho->pipi) by (pipi->rho)*
    FOR_SAMPLE(rs_gevp_corr, s)
    {
        for(int gev = 1; gev < ngev; ++gev)
        {
            rs_gevp_corr[s].col(2*gev) = rs_gevp_corr[s].col(2*gev*ngev);
            rs_gevp_corr[s].col(2*gev + 1) = -rs_gevp_corr[s].col(2*gev*ngev + 1);
        }
    }

    // Generalized eigenvalues
    Sample<Matrix<std::complex<double>>> rs_gevp_gev(params.nboot);
    for(int s = 0; s < params.nboot; ++s)
    {
        rs_gevp_gev[s] = gev(rs_gevp_corr[s], params.t0);
    }
    Matrix<std::complex<double>> gevp_gev = rs_gevp_gev.mean();

    // Effective mass
    out << "computing gev effective masses... ";
    Sample<Matrix<double>> rs_gevp_meff(params.nboot);
    rs_gevp_meff.resizeMatrix(Nt, ngev);
    for(int s = 0; s < params.nboot; ++s) {
        // for(int i = 0; i < t0; ++i) {
        for(int i = 0; i < Nt; ++i) {
            for(int j = 0; j < ngev; ++j)
            {
                rs_gevp_meff[s](i, j) = std::log(fabs(real(rs_gevp_gev[s](i, j)))) / (t0 - i);
            }
        }
        // for(int i = t0; i < Nt; ++i) {
        //     for(int j = 0; j < ngev; ++j)
        //     {
        //         rs_gevp_meff[s](i, j) = std::log(fabs(real(rs_gevp_gev[s](i, (j+1)%2)))) / (t0 - i);
        //     }
        // }
    }

    Matrix<double> gevp_meff = rs_gevp_meff.mean();
    Matrix<double> gevp_meff_err = rs_gevp_meff.variance();
    out << "done\n";

    // Print effective masses
    ofstream of_gevp_meff("gevp_meff.dat");
    for(int i = 0; i < gevp_meff.rows(); ++i) {
        of_gevp_meff << i << ' ';
        for(int j = 0; j < gevp_meff.cols(); ++j)
        {
            of_gevp_meff << gevp_meff(i, j) << ' ' << sqrt(gevp_meff_err(i, j)) << ' ';
        }
        of_gevp_meff << endl;
    }

    // Compute resampled plateaus
    out << "fitting plateaus... \n";
    Sample<Matrix<double>> rs_gevp_plat(rs_gevp_meff.size(), 1, ngev);

    Sample<Matrix<double>> buf;
    for(int j = 0; j < ngev; ++j)
    {
        buf = rs_gevp_meff.col(j);
        auto plat = fitPlateau(buf, params.gevp_range);
        FOR_SAMPLE(rs_gevp_plat, s)
        {
            rs_gevp_plat[s](0, j) = plat[s];
        }
    }
    out << "done\n";

    Matrix<double> gevp_plat = rs_gevp_plat.mean();
    Matrix<double> gevp_plat_err = rs_gevp_plat.variance();

    // Print results
    out << "GEVP analysis : " << endl;
    for(int j = 0; j < ngev; ++j)
    {
        out << "E" << j << "_COM = " << gevp_plat(0, j) << " (" << GeV(gevp_plat(0, j), params.beta) << " GeV)"
            << " +- " << sqrt(gevp_plat_err(0, j)) << " (" << GeV(sqrt(gevp_plat_err(0, j)), params.beta) << " GeV)"
            << endl;
    }
    out << "********************************************************************" << endl << endl;

    out << endl;

    /*** Luscher analysis ***/

    for(int gev1 = 0; gev1 < ngev; ++gev1)
        for(int gev2 = gev1 + 1; gev2 < ngev; ++gev2)
        {
            out << endl
            << "********************************************************************" << endl
            << "Performing Luscher analysis with levels " << gev1 << " and " << gev2 << " :" << endl;

            Sample<double> rs_mrho2(params.nboot);
            Sample<double> rs_g2(params.nboot);
            Sample<double> rs_mrho(params.nboot);
            Sample<double> rs_g(params.nboot);

            FOR_SAMPLE(rs_g2, s)
            {
                double E1 = rs_gevp_plat[s](0, gev1);
                double E2 = rs_gevp_plat[s](0, gev2);
                double E2_1_com = E1*E1;
                double E2_2_com = E2*E2;
                double MPicom = rs_pi_plat[s];
                double DeltaM2_1 = deltaM2_COM(E1, MPicom, params.L);
                double DeltaM2_2 = deltaM2_COM(E2, MPicom, params.L);

                double MRho2 = fabs((E2_1_com*DeltaM2_2 - E2_2_com*DeltaM2_1)/(DeltaM2_2 - DeltaM2_1));
                double g2 = 6.*M_PI*fabs((E2_1_com-E2_2_com)/(DeltaM2_1-DeltaM2_2));

                rs_mrho2[s] = MRho2;
                rs_g2[s] = g2;

                rs_mrho[s] = sqrt(MRho2);
                rs_g[s] = sqrt(g2);
            }

            double mrho2 = rs_mrho2.mean();
            double mrho2_err = rs_mrho2.variance();
            double mrho = rs_mrho.mean();
            double mrho_err = rs_mrho.variance();
            double g2 = rs_g2.mean();
            double g2_err = rs_g2.variance();
            double g = rs_g.mean();
            double g_err = rs_g.variance();

            out << "Coupling parameters : " << endl
            << "g2 = " << g2
            << "  +- " <<  sqrt(g2_err)
            << endl
            << "Mrho2 = " << mrho2
            << "  +- " <<  sqrt(mrho2_err)
            << endl << endl
            << "g = " << g
            << " +- " <<  sqrt(g_err)
            << endl
            << "Mrho = " << mrho << " (" << GeV(mrho, params.beta) << " GeV)"
            << " +- " << mrho_err << " (" << GeV(sqrt(mrho_err), params.beta) << " GeV)"
            << endl;
            out << "********************************************************************" << endl << endl;

    }

 	return 0;
 }