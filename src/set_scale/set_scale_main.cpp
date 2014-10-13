/*
 * set_scale_main.cpp
 *
 *  Created on: Jul 23, 2014
 *      Author: Thibaut Metivet
 */

#include <iostream>
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
    ("L", po::value<unsigned int>(), "space extent of the lattice")
    ("T", po::value<unsigned int>(), "time extent of the lattice")
    ("nboot", po::value<unsigned int>(&nBoot)->default_value(2000), "set number of bootstraps used for statistical resampling analysis")
    ("pi-fit-range", po::value(&pi_fit_range), "set fit range for local ops plateaus")
    ("K-fit-range", po::value(&K_fit_range), "set fit range for gevp ops plateaus")
    ("Omega-fit-range", po::value(&Omega_fit_range), "set fit range for gevp ops plateaus")
    ("Mpi", po::value<double>()->default_value(0.134977), "physical pion mass (GeV)")
    ("MK", po::value<double>()->default_value(0.497614), "physical pion mass (GeV)")
    ("MOmega", po::value<double>()->default_value(1.67245), "physical pion mass (GeV)")
    ;

    po::options_description hidden("Parameters");
    hidden.add_options()
    ("inputfile", po::value<std::string>(), "process a single file");

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

    set_scale_parameters params
    {
        pi_fit_range,
        K_fit_range,
        Omega_fit_range,
        nBoot,
        vm["Mpi"].as<double>(),
        vm["MK"].as<double>(),
        vm["MOmega"].as<double>()
    };

    set_scale_result result = set_scale(vm["inputfile"].as<std::string>(), params);

    cout << "********************************************************************" << endl
         << "Scale setting:" << endl
         << "a = " << result.a << " +- " << result.a_err << " GeV-1" << endl
         << "  = " << result.a * 0.197420751 << " +- " << result.a_err * 0.197420751 << " fm" << endl
         << endl
         << "a^-1 = " << result.a_inv << " +- " << result.a_inv_err << " GeV" << endl
         ;

    return 0;
}

