/*
 * model_modified_luscher.cpp
 *
 *  Created on: Jun 06, 2013
 *      Author: Thibaut Metivet
 */

#include "model_modified_luscher.hpp"
#include <iostream>
#include <gsl/gsl_errno.h>

using namespace LocalLuscher;

double x[]={0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95};
double y[]={-0.282095, -0.395855, -0.484329, -0.545197, -0.57565, -0.572241,-0.530666, -0.445455, -0.309502, -0.113352, 0.155947, 0.516735,0.996053, 1.63577, 2.50485, 3.72621, 5.54399, 8.52506, 14.3652,31.5308};

/*class Q2_Helper
{
private:
    double _L2, _Mpi2, _Mrho2, _g2;
    const LQCDA::FunctionInterpolator& _Q2_Z;

public:
    Q2_Helper(double L2, double Mpi2, double Mrho2, double g2, const LQCDA::FunctionInterpolator& q2_z) :
	_L2(L2), _Mpi2(Mpi2), _Mrho2(Mrho2), _g2(g2), _Q2_Z(q2_z)
	{}

    double operator() (double q2) {
	double tmp = _Mpi2*_L2 + 4.*M_PI*M_PI*q2;
	return ( _Mrho2*_L2 - 4.*tmp ) * sqrt(tmp) * 3. / (2.*M_SQRTPI)
	    - _g2*_Q2_Z(q2);
    }
    };*/

struct q2_helper_params {
    double L2, Mpi2, Mrho2, g2;
    const LQCDA::FunctionInterpolator& Q2_Z;

    q2_helper_params(double l2, double mpi2, double mrho2, double gg2, const LQCDA::FunctionInterpolator& q2_Z) :
	L2(l2), Mpi2(mpi2), Mrho2(mrho2), g2(gg2), Q2_Z(q2_Z) {}
};

double q2_helper(double q2, void* params) {
    q2_helper_params* p = (q2_helper_params*) params;
    double L2 = p->L2;
    double Mpi2 = p->Mpi2;
    double Mrho2 = p->Mrho2;
    double g2 = p->g2;
    const LQCDA::FunctionInterpolator& Q2_Z = p->Q2_Z;

    double tmp = Mpi2*L2 + 4.*M_PI*M_PI*q2;
    return ( Mrho2*L2 - 4.*tmp ) * sqrt(tmp) * 3. / (2.*M_SQRTPI)
	- g2*Q2_Z(q2);
}

double FindRoot(gsl_function* fun, double x_min, double x_max, double epsrel =1.e-4) {
    double r = 0.;
int status;
int iter = 0, max_iter = 100;
double x_lo = x_min;
double x_hi = x_max;

gsl_root_fsolver * Solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
gsl_root_fsolver_set(Solver, fun, x_lo, x_hi);

do {
iter++;
status = gsl_root_fsolver_iterate(Solver);
r = gsl_root_fsolver_root(Solver);

x_lo = gsl_root_fsolver_x_lower(Solver);
x_hi = gsl_root_fsolver_x_upper(Solver);
status = gsl_root_test_interval(x_lo, x_hi, 0., epsrel);

} while(status == GSL_CONTINUE && iter < max_iter);

gsl_root_fsolver_free(Solver);

return r;
}

LocalLuscher::ModifiedLuscherModel::ModifiedLuscherModel() :
    _Q2_Z(x, y,20)
{}

std::vector<double> LocalLuscher::ModifiedLuscherModel::operator() (const std::vector<double>& x, const std::vector<double>& params) const
{
    double L = x[0];
    double Mpi2 = x[1];

    double Mrho = params[0] * Mpi2 + params[1]; // Mrho = a * Mpi2 + b (affine)
    double Mrho2 = Mrho * Mrho;
    double g = params[2];
	
    q2_helper_params* p = new q2_helper_params(L*L, Mpi2, Mrho2, g*g, _Q2_Z);

    double x_lo = 0.;
    double x_hi = 0.95;
    
    if(q2_helper(x_lo, p) * q2_helper(x_hi, p) > 0) {
	return std::vector<double>(1, 1.e4);
    }
	
    gsl_function q2_fun;
    q2_fun.function = &q2_helper;   
    q2_fun.params = p;

    double q2 = FindRoot(&q2_fun, x_lo, x_hi);

    double E = 2. * sqrt( 4.*M_PI*M_PI*q2/(L*L) + Mpi2 );

    delete p;
	    
    return std::vector<double>(1, E);
}
