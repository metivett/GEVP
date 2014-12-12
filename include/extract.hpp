/*
 * extract.hpp
 *
 *  Created on: Sep 15, 2014
 *      Author: Thibaut Metivet
 */

#ifndef EXTRACT_HPP
#define EXTRACT_HPP

#include "LQCDA.hpp"
#include "utils.hpp"

#include "plateau.hpp"
#include "hadron.hpp"

struct extract_parameters
{
    // Particle
    Hadron particle;
    // Smearing
    std::string smearing;
    // number of bootstrap samples
    unsigned int nboot;
    // Beta
    double beta;
    // Delay plateau fit range
    bool delay;
    // Fold correlator at T/2
    bool fold_correlator;
    // Save result
    bool save;
    std::string save_file;
};

inline std::ostream& operator<<(std::ostream& out, const extract_parameters& p)
{
    return out;
}

Hadron hadron_from_string(const std::string &s);

// int extract(const std::string& manfile, const extract_parameters& params);
LQCDA::Sample<LQCDA::Matrix<double>> extract_corr(const std::string& manfile, const extract_parameters& params);
LQCDA::Sample<LQCDA::Matrix<double>> extract_meff(const std::string& manfile, const extract_parameters& params, MeffType mefftype);
LQCDA::Sample<double> extract_mass(const std::string& manfile, const extract_parameters& params, MeffType mefftype, const fit_range& range);
LQCDA::Sample<LQCDA::Matrix<double>> extract_mass_varpro(const std::string &manfile, const extract_parameters &params, MeffType mefftype, unsigned int nExp, const fit_range &range);

#endif // EXTRACT_HPP
