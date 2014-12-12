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

 #include "utils.hpp"


 enum analysis_frame
 {
 	COM = 0,
 	MV = 1
 };

 struct analysis_parameters
 {
 	// lattice dimensions
 	unsigned int L, T;
 	// beta
 	double beta;
    // fold correlators
    bool fold_corr;
 	// analysis frame
 	analysis_frame frame;
 	// GEVP t0
 	unsigned int t0;
 	// number of bootstrap samples
 	unsigned int nboot;
    // pion correlators alternative manfile
    std::string pi_manfile;
 	// local plateau fit range
 	fit_range local_range;
 	// gevp plateau fit range
 	fit_range gevp_range;
    // sin2d vs (E,mpi) fit range
    fit_range sin2d_range;
 };

 std::ostream& operator<<(std::ostream& os, const analysis_parameters& p);

 int analyze(const std::string& manfile, const analysis_parameters& params);

#endif // ANALYZE_HPP