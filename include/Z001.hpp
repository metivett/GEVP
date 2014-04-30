/*
 * Z001.hpp
 *
 *  Created on: Jul 05, 2013
 *      Author: Thibaut Metivet
 */

#ifndef Z001_HPP_
#define Z001_HPP_

 #include <gsl/gsl_math.h>

 namespace Z001 {

 	int n2degen(unsigned int n2);
 	double z001q2(double q2);
 	double z001q2(double q2, double gamma, int direction [3]);

 }

#endif	// Z001_HPP_
