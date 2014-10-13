/*
 * gevp.hpp
 *
 *  Created on: Apr 24, 2014
 *      Author: Thibaut Metivet
 */

#ifndef GEVP_HPP
#define GEVP_HPP

 #include "LQCDA.hpp"
 #include "utils.hpp"

 #include <gsl/gsl_math.h>
 #include <complex>

 LQCDA::Matrix<std::complex<double>> gev(const LQCDA::Matrix<double>& mat, unsigned int t0);
 LQCDA::Matrix<double> gev_svd(const LQCDA::Matrix<double>& mat, unsigned int t0);

 double deltaM2_COM(double W, double MPi, double L);

 // Luscher's cotg(delta)
 double cotd_COM(double W, double MPi, double L);

 namespace Models
 {
 	class Sin2d_vs_E_mpi
 	: public LQCDA::ParametrizedScalarFunction<double>
 	{
 	public:
 		Sin2d_vs_E_mpi()
 		: LQCDA::ParametrizedScalarFunction<double>(2, 2)
 		{}

 		virtual double operator()(const double* x, const double* p) const override
 		{
 			double E = x[0];
 			double E2 = E * E;
 			double mpi2 = x[1] * x[1];
 			double g2 = p[0];
 			double mrho2 = p[1];

 			double k = sqrt(E2/4. - mpi2);
 			double k3 = k*k*k;

 			double tmp = 6.*M_PI/g2 * E/k3 *(mrho2 - E2);
 			return 1. / ( 1. + (tmp*tmp) );
 		}
 	};
 }

 LQCDA::Sample<LQCDA::Matrix<double>> fitSin2d(
 	const LQCDA::Matrix<double>& gevp_plat, 
 	const LQCDA::Sample<LQCDA::Matrix<double>>& rs_sin2d, 
 	const LQCDA::Sample<double>& rs_pi_plat,
    const fit_range& range = fit_range());

#endif // GEVP_HPP