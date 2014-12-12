/*
 * hadron.cpp
 *
 *  Created on: Oct 28, 2014
 *      Author: Thibaut Metivet
 */

#include "hadron.hpp"

std::map<Hadron, hadron_property_sheet> hadron_properties =
{
{Hadron::pion, {138.04, +1, false, false}},
{Hadron::kaon, {494.99, -1, false}},
{Hadron::rho, {775.5, -1, false, false}},

{Hadron::nucleon, {938.92, +1, true, false}},
{Hadron::sigma, {1193.15, +1, true, true}},
{Hadron::xi, {1318.3, +1, true, true}},

{Hadron::omega, {1672.45, +1, true, true}},
{Hadron::xi_star, {1533.4, +1, true, true}},
{Hadron::sigma_star, {1384.6, +1, true, true}},
{Hadron::delta, {1232., +1, true, false}}
};