/*
 * main.cpp
 *
 *  Created on: Apr 16, 2014
 *      Author: Thibaut Metivet
 */

 #include <iostream>
 #include <boost/program_options.hpp>
 #include <string>
 #include <sstream>
 #include <vector>

 #include "LQCDA.hpp"
 #include "analyze.hpp"

 using namespace std;
 namespace po = boost::program_options;

 // void validate(boost::any& v, const vector<string>& values, fit_range*, int)
 // {
 // 	fit_range range;
 // 	if(values.size() != 2)
 // 	{
 // 		throw po::validation_error(
 //            "Invalid fit range specification");
 // 	}
 // 	range.tmin = utils::strTo<int>(values.at(0));
 // 	range.tmax = utils::strTo<int>(values.at(1));
 // }

 std::istream& operator>>(std::istream& src, fit_range& r)
 {
 	string buf, it;
 	src >> buf;
 	stringstream ss(buf);
 	vector<int> range;
 	while(getline(ss, it, ','))
 	{
 		range.push_back(utils::strTo<int>(it));
 	}
 	r.tmin = range.at(0);
 	r.tmax = range.at(1);
 	return src;
 }

 int main(int argc, char **argv) {
 	unsigned int nBootstraps;

 	po::positional_options_description p;
 	p.add("inputfile", -1);

 	po::options_description generic("Generic options");
 	generic.add_options()
 	("help", "print help message");
	// ("inputlist", po::value<std::string>(), "process list of input files")

	fit_range local_range, gevp_range, sin2d_range;	

 	po::options_description parameters("Parameters");
 	parameters.add_options()
 	("t0", po::value<unsigned int>()->required(), "set the t0 reference time")
 	("L", po::value<unsigned int>()->required(), "space extent of the lattice")
 	("T", po::value<unsigned int>()->required(), "time extent of the lattice")
 	("frame", po::value<unsigned int>()->default_value(0), "set type of analysis (COM or MV)")
	// ("t1", po::value<int>(), "if t1 > 0, apply t1-trick exponential correction")
 	("beta", po::value<double>()->default_value(0.), "provide optional beta value to convert lattices units into physical units")
    ("fold-corr", po::value<bool>()->default_value(false), "fold correlators")
 	("nboot", po::value<unsigned int>(&nBootstraps)->default_value(2000), "set number of bootstraps used for statistical resampling analysis")
    ("pi-manfile", po::value<std::string>()->default_value(""), "set alternative manfile for pion correlators")
 	("loc-fit-range", po::value(&local_range), "set fit range for local ops plateaus")
 	("gevp-fit-range", po::value(&gevp_range), "set fit range for gevp ops plateaus")
    ("sin2d-fit-range", po::value(&sin2d_range), "set fit range for sin2d vs (E,mpi) fit")
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

 	if (vm.count("help")) {
 		cout << visible << "\n";
 		return 0;
 	}
 	po::notify(vm);

 	analysis_parameters params {
 		vm["L"].as<unsigned int>(),
 		vm["T"].as<unsigned int>(),
 		vm["beta"].as<double>(),
        vm["fold-corr"].as<bool>(),
 		analysis_frame(vm["frame"].as<unsigned int>()),
 		vm["t0"].as<unsigned int>(),
 		vm["nboot"].as<unsigned int>(),
        vm["pi-manfile"].as<std::string>(),
 		local_range,
 		gevp_range,
        sin2d_range
 	};
 	
 	return analyze(vm["inputfile"].as<std::string>(), params);
 }