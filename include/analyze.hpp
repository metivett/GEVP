/*
 * analyze.hpp
 *
 *  Created on: Apr 16, 2014
 *      Author: Thibaut Metivet
 */

#ifndef ANALYZE_HPP
#define ANALYZE_HPP

 #include <string>
 #include <ostream>
 #include <vector>


 enum analysis_frame
 {
 	COM = 0,
 	MV = 1
 };

 struct fit_range
 {
 	int tmin{0}, tmax{0};
 };

 struct analysis_parameters
 {
 	// lattice dimensions
 	unsigned int L, T;
 	// beta
 	double beta;
 	// analysis frame
 	analysis_frame frame;
 	// GEVP t0
 	unsigned int t0;
 	// number of bootstrap samples
 	unsigned int nboot;
 	// local plateau fit range
 	fit_range local_range;
 	// gevp plateau fit range
 	fit_range gevp_range;
 };

 std::ostream& operator<<(std::ostream& os, const analysis_parameters& p);

 int analyze(const std::string& manfile, const analysis_parameters& params);

#endif // ANALYZE_HPP