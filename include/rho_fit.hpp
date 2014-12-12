/*
 * rho_fit.hpp
 *
 *  Created on: Dec 12, 2014
 *      Author: Thibaut Metivet
 */

#ifndef RHO_FIT_HPP
#define RHO_FIT_HPP

#include "LQCDA.hpp"

enum class FunctionalForm
{
    NONE,
    POLY1,
    POLY2
};

struct rho_fit_parameters
{
    // manfile for data from single correlators
    std::string single_manfile;
    // manfile for data from GEVP correlators
    std::string gevp_manfile;

    // functional form for a dependence
    FunctionalForm a_functional;
    // functional form for L dependence
    FunctionalForm L_functional;
    // functional form for mpi dependence
    FunctionalForm mpi_functional;
    // functional form for mK dependence
    FunctionalForm mK_functional;
};

struct rho_fit_result
{

};

#endif // RHO_FIT_HPP
