/*
 * GEVP.hpp
 *
 *  Created on: Jun 14, 2013
 *      Author: Thibaut Metivet
 */

#ifndef GEVP_HPP_
#define GEVP_HPP_

#include "GEVPHeaders.hpp"
#include "GEVPPlateau.hpp"

std::vector<std::vector<Eigen::Matrix2cd> > ResampleGEVPCorrelators(const std::vector<std::vector<Eigen::Matrix2cd> >& sample, unsigned int Nboot, LQCDA::RandGen::rg_state state);

std::vector<GEVPPlateau> ComputeGEV(const std::vector<Eigen::Matrix2cd>& Cts, int t0);

std::vector<GEVPPlateau> ComputeGEVPPlateau(const std::vector<std::vector<Eigen::Matrix2cd> >& rsplat, int t0, int tf);

std::vector<Eigen::Matrix2cd> ReadGEVPCorrelator(const GEVPCorrelatorHeader& h, const std::string& file);

std::vector<std::vector<Eigen::Matrix2cd> > ReadGEVPSample(const GEVPCorrelatorHeader& h, const std::string& list);

std::vector<std::vector<double> > ResampleLocalCorrelators(const std::vector<std::vector<double> >& sample, unsigned int Nboot, LQCDA::RandGen::rg_state state);

std::vector<double> ReadLocalCorrelator(const GEVPCorrelatorHeader& h, const std::string& file);

std::vector<std::vector<double> > ReadLocalSample(const GEVPCorrelatorHeader& h, const std::string& list);

std::vector<double> ComputeLocalPlateau(const std::vector<std::vector<double> >& rsplat, int tsize);

double DeltaM2_COM(double W, double MPi, double L);
double DeltaM2_MV(double W, double MPi, double L, int d [3]);

#endif	// GEVP_HPP_
