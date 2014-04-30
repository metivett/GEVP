/*
 * GEVPFitModels.hpp
 *
 *  Created on: Jul 01, 2013
 *      Author: Thibaut Metivet
 */

#ifndef GEVP_FIT_MODELS_HPP_
#define GEVP_FIT_MODELS_HPP_

#include "LQCDA.hpp"

    class PlateauFitModel : public LQCDA::FitModel<double, double>
    {
    public:
	double operator() (const double& x, const std::vector<double>& params) const {
	    assert(params.size() == 3);
	    double m = params[0];
	    double E = params[1];
	    double A = params[2];
	    return A*std::exp(- E * x) + m;
	}

	unsigned int NbOfParameters() const { return 3; }
	unsigned int XDim() const { return 1; }
	unsigned int YDim() const { return 1; }
    };

    class ConstantPlateauFitModel : public LQCDA::FitModel<double, double>
    {
    public:
	double operator() (const double& x, const std::vector<double>& params) const {
	    assert(params.size() == 1);
	    double m = params[0];
	    return m;
	}

	unsigned int NbOfParameters() const { return 1; }
	unsigned int XDim() const { return 1; }
	unsigned int YDim() const { return 1; }
    };


#endif	// GEVP_FIT_MODELS_HPP_
