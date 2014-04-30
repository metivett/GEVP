/*
 * model_luscher.hpp
 *
 *  Created on: Apr 18, 2013
 *      Author: Thibaut Metivet
 */

#ifndef MODEL_LUSCHER_HPP_
#define MODEL_LUSCHER_HPP_

#include "LocalLuscher.hpp"
#include <gsl/gsl_math.h>
#include <cmath>

namespace LocalLuscher {
    
    // Function to fit the Luscher formula
    // x = (L,Mpi^2, E^2)
    // y = 0
    // Params = Mrho(a,b), g
    class LuscherModel : public Model
    {
    private:
	static const int N2MAX = 10001;
	static const double EPSILON;
	
	static int _n2degen[N2MAX];

	static int degen(size_t n2);
	
	static double z001q2(double q2);
    public:
	LuscherModel() { build_n2degen(); }

	static void build_n2degen();
	static int n2degen(size_t n2) {
	    if(n2 < N2MAX)
		return _n2degen[n2];
	    else
		return degen(n2);
	}

	virtual unsigned int XDim() const { return 3; }
	virtual unsigned int YDim() const { return 1; }
	virtual unsigned int NbOfParameters() const { return 3; }

	std::vector<double> operator() (const std::vector<double>& x, const std::vector<double>& params) const {
	    double L = x[0];
	    double Mpi2 = x[1];
	    double E = x[2];
	    double E2 = E * E;
	    
	    double k2 = E2/(double)4 - Mpi2;
	    double q2 = 0.25 * k2 * L * L / (M_PI * M_PI);
	    double Dm2 = - 2.0 * k2 * z001q2(q2) / (E * L * M_SQRTPI);
	    
	    double Mrho = params[0] * Mpi2 + params[1]; // Mrho = a * Mpi2 + b (affine)
	    double Mrho2 = Mrho * Mrho;
	    double g = params[2];
	    
	    std::vector<double> result(1, Mrho2 - E2 + g * g * Dm2 / (6 * M_PI));
	    return result;
	}
    };
    
} // namespace LocalLuscher

#endif	// MODEL_LUSCHER_HPP_
