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
namespace po = boost::program_options;

MeffType meff_type_from_str(const string& s)
{
    if(s=="COSH")
        return MeffType::COSH;
    else if(s=="SINH")
        return MeffType::SINH;
    else if(s=="LOG")
        return MeffType::LOG;
    else
        throw std::runtime_error("bad MeffType");
}

int main(int argc, char **argv)
{
    po::options_description generic("Generic options");
    generic.add_options()
    ("help", "print help message");

    po::positional_options_description p;
    p.add("manfile", 1);
    p.add("op", 1);

    unsigned int nBoot;

    po::options_description parameters("Parameters");
    parameters.add_options()
    ("meff", po::value<std::string>(), "extract the effective mass with given method")
    ("save-corr", po::value<std::string>(), "to save the correlator")
    ("save-meff", po::value<std::string>(), "to save the effective mass")
    ("nboot", po::value<unsigned int>(&nBoot)->default_value(2000), "set number of bootstraps used for statistical resampling analysis")
    ;

    po::options_description hidden("Parameters");
    hidden.add_options()
    ("manfile", po::value<std::string>(), "manifest file")
    ("op", po::value<std::string>(), "interpolator to extract");

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

    bool extract_meff = vm.count("meff");
    bool save_corr = vm.count("save-corr");
    bool save_meff = vm.count("save-meff");

    extract_parameters params
    {
        vm["op"].as<string>(),
        vm["nboot"].as<unsigned int>(),
        extract_meff,
        meff_type_from_str(vm["meff"].as<string>()),
        save_corr,
        save_corr? vm["save-corr"].as<string>(): "",
        save_meff,
        save_meff? vm["save-meff"].as<string>(): ""
    };

    return extract(vm["manfile"].as<string>(), params);

    return 0;
}