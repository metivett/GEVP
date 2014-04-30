/*
 * model_modified_luscher.hpp
 *
 *  Created on: Jun 06, 2013
 *      Author: Thibaut Metivet
 */

#ifndef MODEL_MODIFIED_LUSCHER_HPP_
#define MODEL_MODIFIED_LUSCHER_HPP_

#include "LocalLuscher.hpp"
#include "model_luscher.hpp"

namespace LocalLuscher {

// Function to fit the Luscher formula (modified version)
// x = (L,Mpi^2,a)
// y = E
// Params = Mrho(a,b), g
    class ModifiedLuscherModel : public Model
    {
    private:
	LQCDA::FunctionInterpolator _Q2_Z;	// Interpolates q2 * Z(q2)

    public:
	ModifiedLuscherModel();

	virtual unsigned int XDim() const { return 2; }
	virtual unsigned int YDim() const { return 1; }
	virtual unsigned int NbOfParameters() const { return 3; }

	std::vector<double> operator() (const std::vector<double>& x, const std::vector<double>& params) const;
    };

} // namespace LocalLuscher

#endif	// MODEL_MODIFIED_LUSCHER_HPP_
