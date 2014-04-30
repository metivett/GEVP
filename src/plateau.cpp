/*
 * plateau.cpp
 *
 *  Created on: Apr 22, 2014
 *      Author: Thibaut Metivet
 */

#ifndef PLATEAU_CPP
#define PLATEAU_CPP

 #include "plateau.hpp"
 #include "gsl/gsl_cdf.h"

 using namespace std;
 using namespace LQCDA;

 // compute p-value of linear model parameter using t-statistics
 // p is the parameter value
 // p0 the null hypothesis value
 // se the standard error on p
 // ndof is the nb of degrees of freedom
 double pvalue(double p, double p0, double se, unsigned int ndof)
 {
 	double t0 = (p - p0) / se;
 	double p_t_le_t0 = gsl_cdf_tdist_P(t0, ndof);
 	double p_val = 2. * (1 - p_t_le_t0);

 	return p_val;
 }

 struct plat_range
 {
 	int start;
 	int len;
 };

 plat_range findPlateau(const XYDataInterface<double>& data)
 {
 	const double p_val_threshold = 0.2;

 	plat_range res{0, 0};

 	index_t start = 0;
 	index_t len = 3;

 	Models::Line * line = new Models::Line();
 	Models::Constant * constant = new Models::Constant();
 	vector<double> p_init = {0., 0.};

 	Chi2Fit<double, MIN::MIGRAD> Fit(data);
 	Fit.options.verbosity = SILENT;
 	
 	index_t best_start = start;
 	double p_val_min = 1;
 	for(; start < data.nPoints() - len; ++start)
 	{
 		Fit.fitAllPoints(false);
 		Fit.fitPointRange(start, start + len - 1, true);
 		auto fit = Fit.fit(*constant, p_init);

 		double chi2_dof = fit.cost() / fit.nDOF();
 		// cout << "chi2_dof = " << chi2_dof << endl;
 		double x_mean = mean(data.x({start, start + len}, 0));
 		double x_var = variance(data.x({start, start + len}, 0));
 		// cout << "x_var = " << x_var << endl;
 		double p0_p_se = sqrt(chi2_dof*(1./len + x_mean*x_mean/x_var));
 		// double p1_p_se = fit.err(0);
 		double p0_p_val = pvalue(fabs(fit.p(0)), 0., p0_p_se, fit.nDOF());

 		cout << "start = " << start << endl;
 		cout << "p0_p_se = " << p0_p_se << endl;
 		cout << "p_val = " << p0_p_val << endl;

 		if(p0_p_val < p_val_min)
 		{
 			p_val_min = p0_p_val;
 			best_start = start;
 		}
 	}

 	res.start = best_start;

 	len = data.nPoints() - best_start;
 	index_t best_len = len;
 	for(; len > 3; --len)
 	{
 		Fit.fitAllPoints(false);
 		Fit.fitPointRange(best_start, best_start + len - 1, true);
 		auto fit = Fit.fit(*line, p_init);

 		double chi2_dof = fit.cost() / fit.nDOF();
 		double x_var = variance(data.x({best_start, best_start + len - 1}, 0));
 		double p1_p_se = sqrt(chi2_dof / x_var);
 		// double p1_p_se = fit.err(1);
 		double p1_p_val = pvalue(fabs(fit.p(1)), 0., p1_p_se, fit.nDOF());
 		// cout << "chi2_dof = " << chi2_dof << endl;
 		cout << "x_var = " << x_var << endl;
 		cout << "p_val = " << p1_p_val << endl;

 		if(p1_p_val > p_val_threshold);
 			best_len = len;
 	}

 	res.len = best_len;

 	return res;
 }

 Sample<double> fitPlateau(
 	const LQCDA::Sample<LQCDA::Matrix<double>>& rs_meff,
 	const fit_range& f_range)
 {
 	int nboot = rs_meff.size();
 	int npts = rs_meff[0].size();

 	Sample<double> result(nboot);

 	Vector<double> tvec(npts);
 	for(int t = 0; t < npts; ++t)
 	{
 		tvec(t) = t;
 	}

 	// fill XYDataSample
 	XYDataSample<double> * meff_t = new XYDataSample<double>(npts, 1, 1, nboot);
 	for(int n = 0; n < nboot; ++n) {
 		meff_t->x({}, 0)[n] << tvec;
 		meff_t->y({}, 0)[n] << rs_meff[n];
 	}
 	meff_t->setCovFromSample();

 	// search plateau
 	plat_range range;
 	if(f_range.tmin <= 0 && f_range.tmax <= 0)
 	{
	 	cout << "searching plateau...\n";
	 	XYData<double> * meff_t_mean = new XYData<double>(npts, 1, 1);
	 	meff_t_mean->x({}, 0) << tvec;
	 	meff_t_mean->y({}, 0) << rs_meff.mean();
	 	meff_t_mean->yyCov(0, 0) = meff_t->yyCov(0, 0);
	 	range = findPlateau(*meff_t_mean);
	 	cout << "plateau found in range (" << range.start << ", "
	 		<< range.start + range.len - 1 << ")\n";
	}
	else
	{
		range = {f_range.tmin, f_range.tmax - f_range.tmin + 1};
		cout << "using plateau range (" << range.start << ", "
	 		<< range.start + range.len - 1 << ")\n";
	}
	

	Models::Constant * model = new Models::Constant();

 	Chi2Fit<double, MIN::MIGRAD> Fit;
 	Fit.options.verbosity = SILENT;

 	double chi2_dof = 0.;

 	vector<double> E0 = {0.};
 	for(int n = 0; n < nboot; ++n) {
		// // Fit E
		Fit.setData(meff_t->getData(n));	
		Fit.fitPointRange(range.start, range.start + range.len - 1);

		auto fit = Fit.fit(*model, E0);
		chi2_dof += fit.cost() / fit.nDOF();

		// // Store results
 		double E = fabs(fit.parameters()[0]);
 		result[n] = E;
 	}
 	chi2_dof /= nboot;
 	cout << "chi2/dof = " << chi2_dof << endl;

 	delete meff_t;
 	delete model;

 	return result;
 }

#endif // PLATEAU_CPP