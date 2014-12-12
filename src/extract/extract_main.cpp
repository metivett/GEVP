/*
 * extract_main.cpp
 *
 *  Created on: Sep 15, 2014
 *      Author: Thibaut Metivet
 */

#include <iostream>
#include <boost/program_options.hpp>

#include "extract.hpp"

using namespace std;
using namespace LQCDA;
namespace po = boost::program_options;

MeffType meff_type_from_str(const string &s)
{
    if (s == "COSH")
        return MeffType::COSH;
    else if (s == "SINH")
        return MeffType::SINH;
    else if (s == "LOG")
        return MeffType::LOG;
    else
        throw std::runtime_error("bad MeffType");
}

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
    p.add("extract_type", 1);
    p.add("manfile", 1);
    p.add("particle", 1);

    unsigned int nBoot;
    fit_range plat_range;

    po::options_description parameters("Parameters");
    parameters.add_options()
    ("meff-type", po::value<std::string>()->default_value("LOG"), "effective mass extraction method")
    ("smearing", po::value<std::string>()->default_value("GG"), "src/snk smearing")
    ("save", po::value<std::string>()->implicit_value("extract.out"), "save the result")
    ("save-chi2", po::value<std::string>()->implicit_value("extractchi2.out"), "save the chi2 per dof")
    ("nboot", po::value<unsigned int>(&nBoot)->default_value(2000), "number of bootstraps used for statistical resampling analysis")
    ("beta", po::value<double>()->default_value(-1.), "beta")
    ("plat-range", po::value(&plat_range), "plateau fit range")
    ("delay", po::value<bool>()->default_value(false), "delay plateau fit range")
    ("varpro", po::value<unsigned int>()->implicit_value(1), "use variable projection method with provided number of exponentials")
    ("fold-corr", po::value<bool>()->default_value(false), "fold correlator at T/2")
    ;

    po::options_description hidden("Parameters");
    hidden.add_options()
    ("extract_type", po::value<std::string>(), "type of extraction (corr, meff or mass)")
    ("manfile", po::value<std::string>(), "manifest file")
    ("particle", po::value<std::string>(), "particle to extract");

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

    string extract_type = vm["extract_type"].as<string>();
    string manfile = vm["manfile"].as<string>();
    string particle = vm["particle"].as<string>();

    string smearing = vm["smearing"].as<string>();
    MeffType meff_type = meff_type_from_str(vm["meff-type"].as<string>());
    
    bool save = vm.count("save");
    string save_file;
    if (save)
        save_file = vm["save"].as<string>();

    unsigned int nboot = vm["nboot"].as<unsigned int>();
    double beta = vm["beta"].as<double>();
    bool delay = vm["delay"].as<bool>();
    bool fold_corr = vm["fold-corr"].as<bool>();

    extract_parameters params;
    params.particle = hadron_from_string(particle);
    params.smearing = smearing;
    params.nboot = nboot;
    params.beta = beta;
    params.delay = delay;
    params.fold_correlator = fold_corr;
    params.save = save;
    params.save_file = save_file;

    // Print info
    cout << endl << endl
         << "#********************************************************************" << endl
         << "#Extracting with:" << endl
         << "#Manifest file " << manfile << endl
         << "#Parameters " << params << endl
         << "#********************************************************************" << endl;

    if (extract_type == "corr")
    {
        auto rs_corr = extract_corr(manfile, params);
        auto corr = rs_corr.mean();
        auto corr_var = rs_corr.variance();
        for (int i = 0; i < corr.rows(); i++)
        {
            cout << i;
            for (int j = 0; j < corr.cols(); j++)
            {
                cout << ' ' << corr(i, j) << ' ' << sqrt(corr_var(i, j));
            }
            cout << endl;
        }
        if (save) // save corr with error
        {
            ofstream ofile(save_file);
            for (int i = 0; i < corr.rows(); i++)
            {
                ofile << i;
                for (int j = 0; j < corr.cols(); j++)
                {
                    ofile << ' ' << corr(i, j) << ' ' << sqrt(corr_var(i, j));
                }
                ofile << endl;
            }
        }
    }
    else if (extract_type == "meff")
    {
        auto rs_meff = extract_meff(manfile, params, meff_type);
        auto meff = rs_meff.mean();
        auto meff_var = rs_meff.variance();
        for (int i = 0; i < meff.rows(); i++)
        {
            cout << i;
            for (int j = 0; j < meff.cols(); j++)
            {
                cout << ' ' << meff(i, j) << ' ' << sqrt(meff_var(i, j));
            }
            cout << endl;
        }
        if (save) // save meff with errors
        {
            ofstream ofile(save_file);
            for (int i = 0; i < meff.rows(); i++)
            {
                ofile << i;
                for (int j = 0; j < meff.cols(); j++)
                {
                    ofile << ' ' << meff(i, j) << ' ' << sqrt(meff_var(i, j));
                }
                ofile << endl;
            }
        }
    }
    else if (extract_type == "mass")
    {
        if (vm.count("varpro"))
        {
            auto rs_a_E = extract_mass_varpro(manfile, params, meff_type, vm["varpro"].as<unsigned int>(), plat_range);
            auto a_E = rs_a_E.mean();
            auto a_E_var = rs_a_E.variance();
            for(unsigned int j=0; j<a_E.cols(); ++j)
            {
                cout << "a" << j << " = " << a_E(0, j) << " +- " << sqrt(a_E_var(0, j)) << endl;
                cout << "E" << j << " = " << a_E(1, j) << " +- " << sqrt(a_E_var(1, j)) << endl;
                cout << endl;
            }
            if (save) // save resampled mass
            {
                ofstream ofile(save_file);
                FOR_SAMPLE(rs_a_E, s)
                {
                    for(unsigned int j=0; j<rs_a_E.cols(); ++j)
                    {
                        ofile << rs_a_E[s](0, j) << ' ' << rs_a_E[s](1, j) << ' ';
                    }
                    ofile << endl;
                }
            }
        }
        else
        {
            auto rs_mass = extract_mass(manfile, params, meff_type, plat_range);
            auto mass = rs_mass.mean();
            auto mass_var = rs_mass.variance();
            cout << mass << ' ' << sqrt(mass_var);
            cout << endl;
            if (save) // save resampled mass
            {
                ofstream ofile(save_file);
                FOR_SAMPLE(rs_mass, s)
                {
                    ofile << rs_mass[s] << endl;
                }
            }
        }
    }

    return 0;
}