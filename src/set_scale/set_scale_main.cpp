/*
 * set_scale_main.cpp
 *
 *  Created on: Jul 23, 2014
 *      Author: Thibaut Metivet
 */

#include <iostream>
#include <iomanip>
#include <boost/program_options.hpp>

#include "set_scale.hpp"

using namespace std;
namespace po = boost::program_options;

std::istream &operator>>(std::istream &src, fit_range &r)
{
    string buf, it;
    src >> buf;
    stringstream ss(buf);
    vector<int> range;
    while (getline(ss, it, ','))
    {
        range.push_back(utils::strTo<int>(it));
    }
    r.tmin = range.at(0);
    r.tmax = range.at(1);
    return src;
}

int main(int argc, char **argv)
{
    po::options_description generic("Generic options");
    generic.add_options()
    ("help", "print help message");

    po::positional_options_description p;
    p.add("inputfile", -1);

    fit_range pi_fit_range, K_fit_range, Omega_fit_range;
    unsigned int nBoot;

    po::options_description parameters("Parameters");
    parameters.add_options()
    ("beta", po::value<std::vector<double>>(), "beta value(s)")
    ("multi", "enable setting of multiple scales simultaneously")
    ("nboot", po::value<unsigned int>(&nBoot)->default_value(2000), "set number of bootstraps used for statistical resampling analysis")
    ("pi-fit-range", po::value(&pi_fit_range), "set fit range for local ops plateaus")
    ("K-fit-range", po::value(&K_fit_range), "set fit range for gevp ops plateaus")
    ("Omega-fit-range", po::value(&Omega_fit_range), "set fit range for gevp ops plateaus")
    ("Mpi", po::value<double>()->default_value(0.1348), "physical pion mass (GeV)")
    ("MK", po::value<double>()->default_value(0.4942), "physical kaon mass (GeV)")
    ("MOmega", po::value<double>()->default_value(1.67245), "physical omega mass (GeV)")
    ("save-rs", po::value<std::string>()->implicit_value("rs_a.dat"), "save resampled a values to file")
    ;

    po::options_description hidden("Parameters");
    hidden.add_options()
    ("inputfile", po::value<std::vector<std::string>>(), "process a single file");

    po::options_description visible("Usage");
    visible.add(generic).add(parameters);


    po::options_description cmdline_opts;
    cmdline_opts.add(generic).add(parameters).add(hidden);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(cmdline_opts).positional(p).run(), vm);

    if (vm.count("help"))
    {
        cout << visible << "\n";
        return 0;
    }
    po::notify(vm);

    if (vm.count("multi"))
    {
        assert(vm["beta"].as<vector<double>>().size() == vm["inputfile"].as<vector<string>>().size());

        set_multi_scales_parameters params
        {
            vm["beta"].as<vector<double>>(),
            nBoot,
            vm["Mpi"].as<double>(),
            vm["MK"].as<double>(),
            vm["MOmega"].as<double>()
        };

        LQCDA::Sample<LQCDA::Matrix<double>> rs_a = set_multi_scales(vm["inputfile"].as<vector<string>>(), params);
        LQCDA::Sample<LQCDA::Matrix<double>> rs_a_inv(params.nboot);
        FOR_SAMPLE(rs_a_inv, s)
        {
            rs_a_inv[s] = 1. / rs_a[s].array();
        }
        LQCDA::Matrix<double> a = rs_a.mean();
        LQCDA::Matrix<double> a_err = sqrt(rs_a.variance().array());
        LQCDA::Matrix<double> a_inv = rs_a_inv.mean();
        LQCDA::Matrix<double> a_inv_err = sqrt(rs_a_inv.variance().array());

        cout << "********************************************************************" << endl
             << "Scale setting:" << endl;
        for (int nb = 0; nb < params.betas.size(); nb++)
        {
            cout << "a(" << setw(4) << params.betas[nb] << ") = " << a(nb, 0) << " +- " << a_err(nb, 0) << " GeV-1" << endl
                 << "        = " << a(nb, 0) * 0.197420751 << " +- " << a_err(nb, 0) * 0.197420751 << " fm" << endl
                 << endl
                 << "a^-1(" << params.betas[nb] <<") = " << a_inv(nb, 0) << " +- " << a_inv_err(nb, 0) << " GeV" << endl
                 ;
        }
        cout << "********************************************************************" << endl;

        if (vm.count("save-rs"))
        {
            ofstream of_rs_a(vm["save-rs"].as<string>());
            FOR_SAMPLE(rs_a, s)
            {
                of_rs_a << rs_a[s].transpose() << endl;
            }
        }
    }
    else
    {
        assert(vm["beta"].as<vector<double>>().size() == 1);
        assert(vm["inputfile"].as<vector<string>>().size() == 1);

        set_scale_parameters params
        {
            pi_fit_range,
            K_fit_range,
            Omega_fit_range,
            vm["beta"].as<vector<double>>()[0],
            nBoot,
            vm["Mpi"].as<double>(),
            vm["MK"].as<double>(),
            vm["MOmega"].as<double>()
        };

        LQCDA::Sample<double> rs_a = set_scale(vm["inputfile"].as<vector<string>>()[0], params);
        LQCDA::Sample<double> rs_a_inv(params.nboot);
        FOR_SAMPLE(rs_a_inv, s)
        {
            rs_a_inv[s] = 1. / rs_a[s];
        }

        double a = rs_a.mean();
        double a_err = sqrt(rs_a.variance());
        double a_inv = rs_a_inv.mean();
        double a_inv_err = sqrt(rs_a_inv.variance());

        cout << "********************************************************************" << endl
             << "Scale setting:" << endl
             << "a = " << a << " +- " << a_err << " GeV-1" << endl
             << "  = " << a * 0.197420751 << " +- " << a_err * 0.197420751 << " fm" << endl
             << endl
             << "a^-1 = " << a_inv << " +- " << a_inv_err << " GeV" << endl
             << "********************************************************************" << endl
             ;

        if (vm.count("save-rs"))
        {
            ofstream of_rs_a(vm["save-rs"].as<string>());
            FOR_SAMPLE(rs_a, s)
            {
                of_rs_a << rs_a[s] << endl;
            }
        }
    }

    return 0;
}

