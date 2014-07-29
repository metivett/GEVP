/*
 * utils.hpp
 *
 *  Created on: Apr 11, 2014
 *      Author: Thibaut Metivet
 */

#ifndef UTILS_HPP
#define UTILS_HPP

 #include <memory>
 #include <string>

 #include "LQCDA.hpp"

 struct fit_range
 {
 	int tmin{0}, tmax{0};
 };

 // convert lattice units to physical units multiplying by a^-1
 double GeV(double x, double beta);

 // read samples
 void readSamples(
 	LQCDA::Sample<LQCDA::Matrix<double>> * sample, 
 	const std::string& manfile,
 	const std::string& header);

 // resample samples
 LQCDA::Sample<LQCDA::Matrix<double>> resample(
 	const LQCDA::Sample<LQCDA::Matrix<double>>& sample, 
 	unsigned int nboot,
 	LQCDA::RandGen::rg_state state);
 

#endif // UTILS_HPP