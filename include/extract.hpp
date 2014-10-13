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

struct extract_parameters
{
    // Correlator header
    std::string corr_header;
    // number of bootstrap samples
    unsigned int nboot;
    // Effective mass
    bool extract_meff;
    MeffType meff_type;
    // Save correlator
    bool save_corr;
    std::string corr_file;
    // Save effective mass
    bool save_meff;
    std::string meff_file;
};

int extract(const std::string& manfile, const extract_parameters& params);

#endif // EXTRACT_HPP
