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

 // fit plateau
 LQCDA::Sample<double> fitPlateau(
 	const LQCDA::Sample<LQCDA::Matrix<double>>& rs_meff,
 	const fit_range& range = fit_range());


#endif // PLATEAU_HPP