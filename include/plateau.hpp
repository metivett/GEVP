/*
 * plateau.hpp
 *
 *  Created on: Apr 22, 2014
 *      Author: Thibaut Metivet
 */

#ifndef PLATEAU_HPP
#define PLATEAU_HPP

 #include "LQCDA.hpp"

 #include "analyze.hpp"

 namespace Models {

 	class Constant
 	: public LQCDA::ParametrizedScalarFunction<double>
 	{
 	public:
 		Constant()
 		: LQCDA::ParametrizedScalarFunction<double>(1, 1)
 		{}

 		virtual double operator()(const double* x, const double* p) const override
 		{
 			return (*p);
 		}
 	};

 	class Line
 	: public LQCDA::ParametrizedScalarFunction<double>
 	{
 	public:
 		Line()
 		: LQCDA::ParametrizedScalarFunction<double>(1, 2)
 		{}

 		virtual double operator()(const double* x, const double* p) const override
 		{
 			return p[0] + p[1] * (*x);
 		}
 	};

 }

 enum class MeffType
 {
 	LOG,
 	COSH,
    SINH
 };

 // local effective mass
 LQCDA::Sample<LQCDA::Matrix<double>> localMeff(
 	MeffType type,
 	const LQCDA::Sample<LQCDA::Matrix<double>>& rs_corr, int Nt);
 //gevp effective mass
 LQCDA::Sample<LQCDA::Matrix<double>> gevpMeff(
 	MeffType type,
 	const LQCDA::Sample<LQCDA::Matrix<std::complex<double>>>& rs_gev, int Nt, int t0);

 // fit plateau
 LQCDA::Sample<double> fitPlateau(
 	const LQCDA::Sample<LQCDA::Matrix<double>>& rs_meff,
 	const fit_range& range = fit_range());

 LQCDA::Sample<double> fitPlateau(
 	const LQCDA::Sample<LQCDA::Matrix<double>>& rs_meff,
 	const fit_range& f_range,
 	LQCDA::Vector<bool>& is_valid);

 // get plateau (from correlator)
 LQCDA::Sample<double> getPlateau(
    const LQCDA::Sample<LQCDA::Matrix<double>>& rs_corr,
    MeffType type,
    const fit_range& range = fit_range());


#endif // PLATEAU_HPP