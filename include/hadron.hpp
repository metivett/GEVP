/*
 * hadron.hpp
 *
 *  Created on: Oct 28, 2014
 *      Author: Thibaut Metivet
 */

#ifndef HADRON_HPP
#define HADRON_HPP

#include <map>
#include <string>

enum class Hadron { pion, kaon, rho, nucleon, sigma, xi, omega, xi_star, sigma_star, delta };

struct hadron_property_sheet
{
    double mass;
    int parity;
    bool is_baryon;
    bool is_strange;
};

extern std::map<Hadron, hadron_property_sheet> hadron_properties;

#endif // HADRON_HPP
