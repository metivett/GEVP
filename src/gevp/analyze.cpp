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
#include <boost/filesystem.hpp>

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
namespace fs = boost::filesystem;

std::ostream &operator<<(std::ostream &os, const analysis_parameters &p)
{
    os
            << "L = " << p.L << endl
            << "T = " << p.T << endl
            << "beta = " << p.beta << endl
            << "frame = " << p.frame << endl
            << "t0 = " << p.t0 << endl
            << "nboot = " << p.nboot << endl
            << "local fit range = {" << p.local_range.tmin << ", " << p.local_range.tmax << "}" << endl
            << "gevp fit range = {" << p.gevp_range.tmin << ", " << p.gevp_range.tmax << "}" << endl
            << "sin2d fit range = {" << p.sin2d_range.tmin << ", " << p.sin2d_range.tmax << "}" << endl;
    return os;
}

int analyze(const std::string &manfile, const analysis_parameters &params)
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
    std::string pi_manfile;
    if (params.pi_manfile != "")
        pi_manfile = params.pi_manfile;
    else
        pi_manfile = manfile;

    readSamples(&pi_corr_sample, pi_manfile, pi_h);
    out << "done\n";

    Matrix<double> pi_t = pi_corr_sample[0].col(0);
    pi_corr_sample = pi_corr_sample.block(0, 1, pi_corr_sample.rows(), pi_corr_sample.cols() - 1);

    // Resample PP correlators
    out << "resampling... ";
    Sample<Matrix<double>> rs_pi_corr = resample(pi_corr_sample, params.nboot, state);
    out << "done\n";

    Vector<double> pi_corr = rs_pi_corr.mean().col(0);
    Vector<double> pi_corr_err = rs_pi_corr.variance().col(0);

    // Print correlator
    ofstream of_pi_corr("pi_corr.dat");
    for (int i = 0; i < pi_corr.size(); ++i)
    {
        of_pi_corr << i << ' ' << pi_corr[i] << ' ' << sqrt(pi_corr_err[i]) << endl;
    }

    // Effective mass
    out << "computing effective mass... ";
    // Sample<Matrix<double>> rs_pi_meff(params.nboot);
    // rs_pi_meff.resizeMatrix(Nt, 1);
    // for(int s = 0; s < params.nboot; ++s) {
    //  for(int i = 0; i < Nt; ++i) {
    //      rs_pi_meff[s](i) = (std::log(rs_pi_corr[s](i, 0) / rs_pi_corr[s](i+1, 0)));
    //  }
    // }
    Sample<Matrix<double>> rs_pi_meff = localMeff(MeffType::COSH, rs_pi_corr, params.T);

    Vector<double> pi_meff = rs_pi_meff.mean();
    Vector<double> pi_meff_err = rs_pi_meff.variance();
    out << "done\n";

    // Print effective mass
    ofstream of_pi_meff("pi_meff.dat");
    for (int i = 0; i < pi_meff.size(); ++i)
    {
        of_pi_meff << i << ' ' << pi_meff[i] << ' ' << sqrt(pi_meff_err[i]) << endl;
    }

    // Compute resampled plateaus
    out << "fitting plateau... \n";
    Vector<bool> rs_pi_plat_validity(params.nboot);
    Sample<double> rs_pi_plat = fitPlateau(rs_pi_meff, params.local_range, rs_pi_plat_validity);
    // Sample<double> rs_pi_plat = fitPlateau(rs_pi_meff, params.local_range);
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
    rho_corr_sample = rho_corr_sample.block(0, 1, rho_corr_sample.rows(), rho_corr_sample.cols() - 1);

    // Resample V1V1 correlators
    out << "resampling... ";
    Sample<Matrix<double>> rs_rho_corr = resample(rho_corr_sample, params.nboot, state);
    out << "done\n";

    Vector<double> rho_corr = rs_rho_corr.mean().col(0);
    Vector<double> rho_corr_err = rs_rho_corr.variance().col(0);

    // Print correlator
    ofstream of_rho_corr("rho_corr.dat");
    for (int i = 0; i < rho_corr.size(); ++i)
    {
        of_rho_corr << i << ' ' << rho_corr[i] << ' ' << sqrt(rho_corr_err[i]) << endl;
    }

    // Effective mass
    out << "computing effective mass... ";
    // Sample<Matrix<double>> rs_rho_meff(params.nboot);
    // rs_rho_meff.resizeMatrix(Nt, 1);
    // for(int s = 0; s < params.nboot; ++s) {
    //     for(int i = 0; i < Nt; ++i) {
    //         rs_rho_meff[s](i) = (std::log(rs_rho_corr[s](i, 0) / rs_rho_corr[s](i+1, 0)));
    //     }
    // }
    Sample<Matrix<double>> rs_rho_meff = localMeff(MeffType::LOG, rs_rho_corr, params.T);

    Vector<double> rho_meff = rs_rho_meff.mean();
    Vector<double> rho_meff_err = rs_rho_meff.variance();
    out << "done\n";

    // Print effective mass
    ofstream of_rho_meff("rho_meff.dat");
    for (int i = 0; i < rho_meff.size(); ++i)
    {
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

    out << "initializing... \n";

    // Read ViVi correlators for all configurations
    std::string gevp_h("ViVi");

    out << "reading samples... ";
    Sample<Matrix<double>> gevp_corr_sample;
    readSamples(&gevp_corr_sample, manfile, gevp_h);
    out << "done\n";

    Matrix<double> gevp_t = gevp_corr_sample[0].col(0);
    gevp_corr_sample = gevp_corr_sample.block(0, 1, gevp_corr_sample.rows(), gevp_corr_sample.cols() - 1);

    int ngev = sqrt((gevp_corr_sample[0].cols()) / 2);
    out << "ngev=" << ngev << endl;

    // Resample ViVi correlators
    out << "resampling... ";
    Sample<Matrix<double>> rs_gevp_corr = resample(gevp_corr_sample, params.nboot, state);
    out << "done\n";

    // Diagonal meffs
    out << "computing pipi effective masses... ";
    Matrix<double> pipi_corr, pipi_meff;
    Matrix<double> pipi_corr_err, pipi_meff_err;
    for (int gev = 1; gev < ngev; ++gev)
    {
        // Sample<Matrix<double>> rs_pipi_corr(rs_gevp_corr.size());
        // rs_pipi_corr.resizeMatrix(rs_gevp_corr.rows(), 1);
        // rs_pipi_corr.col(0) = rs_gevp_corr.col(2*gev*(ngev+1));
        Sample<Matrix<double>> rs_pipi_corr;
        rs_pipi_corr = rs_gevp_corr.col(2 * gev * (ngev + 1));
        Sample<Matrix<double>> rs_pipi_meff = localMeff(MeffType::LOG, rs_pipi_corr, params.T);

        pipi_corr = rs_pipi_corr.mean();
        pipi_corr_err = rs_pipi_corr.variance();
        ofstream of_pipi_corr("pipi" + utils::strFrom(gev) + "_corr.dat");
        for (int i = 0; i < pipi_corr.size(); ++i)
        {
            of_pipi_corr << i << ' ' << pipi_corr(i, 0) << ' ' << sqrt(pipi_corr_err(i, 0)) << endl;
        }

        pipi_meff = rs_pipi_meff.mean();
        pipi_meff_err = rs_pipi_meff.variance();
        ofstream of_pipi_meff("pipi" + utils::strFrom(gev) + "_meff.dat");
        for (int i = 0; i < pipi_meff.size(); ++i)
        {
            of_pipi_meff << i << ' ' << pipi_meff(i, 0) << ' ' << sqrt(pipi_meff_err(i, 0)) << endl;
        }
    }
    out << "done\n";

    // Replace (rho->pipi) by (pipi->rho)*
    FOR_SAMPLE(rs_gevp_corr, s)
    {
        for (int gev = 1; gev < ngev; ++gev)
        {
            rs_gevp_corr[s].col(2 * gev) = rs_gevp_corr[s].col(2 * gev * ngev);
            rs_gevp_corr[s].col(2 * gev + 1) = -rs_gevp_corr[s].col(2 * gev * ngev + 1);
        }
    }
    // Fold (rho->rho) and (pipi->pipi) correlators if requested
    if (params.fold_corr)
    {
        FOR_SAMPLE(rs_gevp_corr, s)
        {
            // (rho->rho)
            for (int t = 1; t < params.T / 2; ++t)
            {
                rs_gevp_corr[s].col(0)(t) = (rs_gevp_corr[s].col(0)(t) + rs_gevp_corr[s].col(0)(params.T - t - 1)) / 2.;
                rs_gevp_corr[s].col(1)(t) = (rs_gevp_corr[s].col(1)(t) + rs_gevp_corr[s].col(1)(params.T - t - 1)) / 2.;
            }
            // (pipi->pipi)
            for (int t = 1; t < params.T / 2; ++t)
            {
                rs_gevp_corr[s].col(2 * (ngev * ngev - 1))(t) = (rs_gevp_corr[s].col(2 * (ngev * ngev - 1))(t) + rs_gevp_corr[s].col(2 * (ngev * ngev - 1))(params.T - t - 1)) / 2.;
                rs_gevp_corr[s].col(2 * ngev * ngev - 1)(t) = (rs_gevp_corr[s].col(2 * ngev * ngev - 1)(t) + rs_gevp_corr[s].col(2 * ngev * ngev - 1)(params.T - t - 1)) / 2.;
            }
        }
    }
    out << "initializing done\n";

    auto t0 = params.t0;

    out << endl;
    // out << "--------------------------------------------------------------------" << endl;
    out << "using t0 = " << t0 << endl;

    // Generalized eigenvalues
    Sample<Matrix<std::complex<double>>> rs_gevp_gev(params.nboot);

    for (int s = 0; s < params.nboot; ++s)
    {
        rs_gevp_gev[s] = gev(rs_gevp_corr[s], t0);
    }
    // Matrix<std::complex<double>> gevp_gev = rs_gevp_gev.mean();

    // Effective mass
    out << "computing gev effective masses... ";
    Sample<Matrix<double>> rs_gevp_meff = gevpMeff(MeffType::LOG, rs_gevp_gev, params.T, t0);

    Matrix<double> gevp_meff = rs_gevp_meff.mean();
    Matrix<double> gevp_meff_err = rs_gevp_meff.variance();
    out << "done\n";

    // Print effective masses
    ofstream of_gevp_meff(utils::strFrom("gevp_meff_t0_", t0, ".dat"));
    for (int i = 0; i < gevp_meff.rows(); ++i)
    {
        of_gevp_meff << i << ' ';
        for (int j = 0; j < gevp_meff.cols(); ++j)
        {
            of_gevp_meff << gevp_meff(i, j) << ' ' << sqrt(gevp_meff_err(i, j)) << ' ';
        }
        of_gevp_meff << endl;
    }

    // Compute resampled plateaus
    out << "fitting plateaus... ";
    Sample<Matrix<double>> rs_gevp_plat(rs_gevp_meff.size(), 1, ngev);

    Sample<Matrix<double>> buf;
    for (int j = 0; j < ngev; ++j)
    {
        buf = rs_gevp_meff.col(j);
        auto plat = fitPlateau(buf, params.gevp_range);
        FOR_SAMPLE(rs_gevp_plat, s)
        {
            rs_gevp_plat[s](0, j) = plat[s];
        }

        fs::path chi2_file("chi2.boot");
        fs::path cur_chi2_file(utils::strFrom("chi2_", j, "_t0_", t0, ".boot"));
        if(fs::exists(chi2_file))
        {
            fs::rename(chi2_file, cur_chi2_file);
        }
    }
    out << "done\n";


    Matrix<double> gevp_plat = rs_gevp_plat.mean();
    Matrix<double> gevp_plat_err = rs_gevp_plat.variance();

    // Print results
    out << "GEVP analysis (t0 = " << t0 << "): " << endl;
    for (int j = 0; j < ngev; ++j)
    {
        out << "E" << j << "_COM = " << gevp_plat(0, j) << " (" << GeV(gevp_plat(0, j), params.beta) << " GeV)"
            << " +- " << sqrt(gevp_plat_err(0, j)) << " (" << GeV(sqrt(gevp_plat_err(0, j)), params.beta) << " GeV)"
            << endl;
    }
    out << "********************************************************************" << endl << endl;

    out << endl;

    /*** Luscher analysis ***/

    out << endl
        << "********************************************************************" << endl
        << "Computing cotg(delta) with Luscher formula: " << endl;
    out << "E\tq\tq2\tcotd\tcotd_err\tsin2d\tsin2d_err" << endl;

    ofstream of_gevp_cotd(utils::strFrom("gevp_cotd_t0", t0, ".dat"));

    of_gevp_cotd << "# E\tq\tq2\tcotd\tcotd_err\tsin2d\tsin2d_err" << endl;

    Sample<Matrix<double>> rs_sin2d(params.nboot, ngev, 1);

    for (int gev = 0; gev < ngev; ++gev)
    {
        Sample<double> rs_cotd(params.nboot);

        ofstream of_rs_cotd(utils::strFrom("rs_cotd", gev, "t0_", t0, ".dat"));
        of_rs_cotd << "# E\tcotd\tsin2d" << endl;

        FOR_SAMPLE(rs_cotd, s)
        {
            double E = rs_gevp_plat[s](0, gev);
            double MPicom = rs_pi_plat[s];
            rs_cotd[s] = cotd_COM(E, MPicom, params.L);
            rs_sin2d[s](gev, 0) = 1. / (1. + rs_cotd[s] * rs_cotd[s]);
            of_rs_cotd << E << " " << rs_cotd[s] << " " << rs_sin2d[s](gev, 0) << endl;
        }

        double cotd = rs_cotd.mean();
        double cotd_err = rs_cotd.variance();
        double sin2d = rs_sin2d.mean()(gev, 0);
        double sin2d_var = rs_sin2d.variance()(gev, 0);
        double q = sqrt((gevp_plat(0, gev) * gevp_plat(0, gev) / 4. - pi_plat * pi_plat) * (params.L * params.L) / (4.*M_PI * M_PI));

        // Print cotg(delta)
        out << gevp_plat(0, gev) << " " << q << " " << q *q << " " << cotd << " " << sqrt(cotd_err) << " " << sin2d << " " << sqrt(sin2d_var) << endl;
        of_gevp_cotd << gevp_plat(0, gev) << " " << q << " " << q *q << " " << cotd << " " << sqrt(cotd_err) << " " << sin2d << " " << sqrt(sin2d_var) << endl;
    }

    for (int gev1 = 0; gev1 < ngev; ++gev1)
        for (int gev2 = gev1 + 1; gev2 < ngev; ++gev2)
        {
            out << endl
                << "--------------------------------------------------------------------" << endl
                << "Performing Luscher analysis with levels " << gev1 << " and " << gev2 << " :" << endl;

            Sample<double> rs_mrho2(params.nboot);
            Sample<double> rs_g2(params.nboot);
            Sample<double> rs_mrho(params.nboot);
            Sample<double> rs_g(params.nboot);

            FOR_SAMPLE(rs_g2, s)
            {
                double E1 = rs_gevp_plat[s](0, gev1);
                double E2 = rs_gevp_plat[s](0, gev2);
                double E2_1_com = E1 * E1;
                double E2_2_com = E2 * E2;
                double MPicom = rs_pi_plat[s];
                double DeltaM2_1 = deltaM2_COM(E1, MPicom, params.L);
                double DeltaM2_2 = deltaM2_COM(E2, MPicom, params.L);

                double MRho2 = fabs((E2_1_com * DeltaM2_2 - E2_2_com * DeltaM2_1) / (DeltaM2_2 - DeltaM2_1));
                double g2 = 6.*M_PI * fabs((E2_1_com - E2_2_com) / (DeltaM2_1 - DeltaM2_2));

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
            out << "--------------------------------------------------------------------" << endl << endl;

            ofstream of_rs_g(utils::strFrom("rs_g", gev1, gev2, "t0_", t0, ".dat"));
            ofstream of_rs_mrho(utils::strFrom("rs_mrho", gev1, gev2, "t0_", t0, ".dat"));
            FOR_SAMPLE(rs_g, s)
            {
                of_rs_g << rs_g[s] << endl;
                of_rs_mrho << rs_mrho[s] << endl;
            }
        }

    out << "********************************************************************" << endl;

    out << endl
        << "--------------------------------------------------------------------" << endl
        << "Performing Luscher analysis with sin^2(delta) fit :" << endl;
    Sample<Matrix<double>> rs_g2_mrho2_from_fit = fitSin2d(gevp_plat, rs_sin2d, rs_pi_plat, params.sin2d_range);

    double rel_err_threshold = 1.;
    // We keep only bootstraps with relative fit errors < rel_err_threshold
    std::vector<unsigned int> usable_boot;
    usable_boot.reserve(params.nboot);
    FOR_SAMPLE(rs_g2_mrho2_from_fit, s)
    {
        double g2_rel_err = rs_g2_mrho2_from_fit[s](1, 0) / rs_g2_mrho2_from_fit[s](0, 0);
        double mrho2_rel_err = rs_g2_mrho2_from_fit[s](1, 1) / rs_g2_mrho2_from_fit[s](0, 1);
        if (g2_rel_err < rel_err_threshold && mrho2_rel_err < rel_err_threshold)
            usable_boot.push_back(s);
    }

    // Sample<double> rs_g_from_fit(usable_boot.size());
    // Sample<double> rs_mrho_from_fit(usable_boot.size());
    // FOR_SAMPLE(rs_g_from_fit, i)
    // {
    //     rs_g_from_fit[i] = sqrt(rs_g2_mrho2_from_fit[usable_boot[i]](0, 0));
    //     rs_mrho_from_fit[i] = sqrt(rs_g2_mrho2_from_fit[usable_boot[i]](0, 1));
    // }
    Sample<double> rs_g_from_fit(rs_g2_mrho2_from_fit.size());
    Sample<double> rs_mrho_from_fit(rs_g2_mrho2_from_fit.size());
    FOR_SAMPLE(rs_g_from_fit, s)
    {
        rs_g_from_fit[s] = sqrt(rs_g2_mrho2_from_fit[s](0, 0));
        rs_mrho_from_fit[s] = sqrt(rs_g2_mrho2_from_fit[s](0, 1));
    }

    double g_from_fit = rs_g_from_fit.mean();
    double g_err_from_fit = rs_g_from_fit.variance();
    double mrho_from_fit = rs_mrho_from_fit.mean();
    double mrho_err_from_fit = rs_mrho_from_fit.variance();

    out << "Coupling parameters from sin^2(delta) fit : " << endl
        << usable_boot.size() << " bootstraps usable" << endl
        << "g = " << g_from_fit
        << " +- " <<  sqrt(g_err_from_fit)
        << endl
        << "Mrho = " << mrho_from_fit << " (" << GeV(mrho_from_fit, params.beta) << " GeV)"
        << " +- " << sqrt(mrho_err_from_fit) << " (" << GeV(sqrt(mrho_err_from_fit), params.beta) << " GeV)"
        << endl;
    out << "--------------------------------------------------------------------" << endl << endl;

    ofstream of_rs_g_from_fit(utils::strFrom("rs_g_from_sin2dt0_", t0, ".dat"));
    ofstream of_rs_mrho_from_fit(utils::strFrom("rs_mrho_from_sin2dt0_", t0, ".dat"));
    FOR_SAMPLE(rs_g_from_fit, s)
    {
        of_rs_g_from_fit << rs_g_from_fit[s] << endl;
        of_rs_mrho_from_fit << rs_mrho_from_fit[s] << endl;
    }

    out << "********************************************************************" << endl;

    return 0;
}