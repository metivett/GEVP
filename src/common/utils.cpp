/*
 * utils.cpp
 *
 *  Created on: Apr 11, 2014
 *      Author: Thibaut Metivet
 */

 #include "utils.hpp"

 #include <map>
 #include <boost/filesystem.hpp>

 using namespace std;
 using namespace LQCDA;
 namespace fs = boost::filesystem;

 double GeV(double x, double beta)
 {
    // Map to convert lattice units to physical units
 	std::map<double, double> ainv{{3.31, 1.697}, {3.5, 2.131}, {3.61, 2.561}, {3.7, 3.026}, {3.8, 3.662}};

 	return x*ainv[beta];
 }

 void readSamples(
 	Sample<Matrix<double>> * sample, 
 	const std::string& manfile,
 	const std::string& header)
 {
 	std::string filename;
 	AsciiDataFile file;

 	std::vector<Matrix<double>> buffer;
 	buffer.reserve(300);

    // std::cout << manfile << std::endl;
 	std::ifstream ifs(manfile);
 	while(std::getline(ifs, filename))
 	{
        fs::path filepath(filename);
        if(!filepath.is_absolute())
        {
            fs::path manfilepath(manfile);
            filepath = manfilepath.parent_path() / filepath;
            filename = filepath.native();
        }

 		file.open(filename, 'r');
 		buffer.push_back(file.getData(header));
 		file.close();
 	}

 	sample->resize(buffer.size());
 	for(int s = 0; s < buffer.size(); ++s)
 	{
 		(*sample)[s] = buffer[s];
 	}
 }

 Sample<Matrix<double>> resample(
 	const LQCDA::Sample<LQCDA::Matrix<double>>& sample, 
 	unsigned int nboot,
 	LQCDA::RandGen::rg_state state)
 {
 	RandGen rng;
 	rng.setState(state);

 	Sample<Matrix<double>> result(nboot);

 	unsigned int NSample = sample.size();
 	double NSample_inv = 1. / (double)NSample;
 	unsigned int index;

 	unsigned int nRow = sample[0].rows();
 	unsigned int nCol = sample[0].cols();

 	for(int n = 0; n < nboot; ++n) {
 		index = rng.getUniformInt(NSample);
 		// Matrix<double> buf(sample[index].block(0, 1, nRow, nCol-1));
        Matrix<double> buf(sample[index]);
 		for(int i = 1; i < NSample; ++i) {
 			index = rng.getUniformInt(NSample);
 			// buf = buf + sample[index].block(0, 1, nRow, nCol-1);
            buf = buf + sample[index];
 		}
 		buf = buf * NSample_inv;

 		result[n] = buf;
 	}

 	return result;
 }
