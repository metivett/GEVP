/*
 * gevp.hpp
 *
 *  Created on: Apr 24, 2014
 *      Author: Thibaut Metivet
 */

#ifndef GEVP_HPP
#define GEVP_HPP

 #include "LQCDA.hpp"

 #include <complex>

 LQCDA::Matrix<std::complex<double>> gev(const LQCDA::Matrix<double>& mat, unsigned int t0);
 LQCDA::Matrix<double> gev_svd(const LQCDA::Matrix<double>& mat, unsigned int t0);

 double deltaM2_COM(double W, double MPi, double L);

#endif // GEVP_HPP