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

class CoshMeffHelper
    : public ScalarFunction<double>
{
private:
    double _t1, _t2, _C;

public:
    CoshMeffHelper(double t1, double t2, double C)
        : ScalarFunction<double>(1)
        , _t1(t1)
        , _t2(t2)
        , _C(C)
    {}

    virtual double operator()(const double *meff) const override
    {
        // cout << "t1 = " << _t1 << endl;
        // cout << "t2 = " << _t2 << endl;
        // cout << "C = " << _C << endl;
        // cout << "meff = " << *meff << endl;
        // cout << _C - cosh(*meff*_t1)/cosh(*meff*_t2) << endl;
        return _C - cosh(*meff * _t1) / cosh(*meff * _t2);
    }
};
class SinhMeffHelper
    : public ScalarFunction<double>
{
private:
    double _t1, _t2, _C;

public:
    SinhMeffHelper(double t1, double t2, double C)
        : ScalarFunction<double>(1)
        , _t1(t1)
        , _t2(t2)
        , _C(C)
    {}

    virtual double operator()(const double *meff) const override
    {
        // cout << "t1 = " << _t1 << endl;
        // cout << "t2 = " << _t2 << endl;
        // cout << "C = " << _C << endl;
        // cout << "meff = " << *meff << endl;
        // cout << _C - sinh(*meff*_t1)/sinh(*meff*_t2) << endl;
        return _C - sinh(*meff * _t1) / sinh(*meff * _t2);
    }
};

fit_range automatic_plat_range(double beta, unsigned int T, bool is_strange, bool delayed)
{
    int ti;
    fit_range plat_range;

    if (beta == 3.31)
        ti = delayed ? 10 : 9;
    else if (beta == 3.5)
        ti = delayed ? 11 : 10;
    else if (beta == 3.61)
        ti = delayed ? 14 : 12;
    else if (beta == 3.7)
        ti = delayed ? 15 : 13;
    else if (beta == 3.8)
        ti = delayed ? 16 : 14;

    int tf = std::min((int)(2.*ti + 0.5), (int)(T / 2 - 2));

    if (is_strange)
        ti *= 1.3;

    fit_range result;
    result.tmin = ti; result.tmax = tf;
    return result;
}

Sample<Matrix<double>> localMeff(
                        MeffType type,
                        const Sample<Matrix<double>> &rs_corr, int Nt)
{
    int nboot = rs_corr.size();
    Sample<Matrix<double>> rs_meff(nboot);

    if (type == MeffType::COSH)
    {
        Roots::BrentRootFinder<double> solver;

        rs_meff.resizeMatrix(Nt - 2, 1);
        for (int s = 0; s < nboot; ++s)
        {
            for (int i = 0; i < Nt - 2; ++i)
            {
                double tmp = fabs(rs_corr[s](i, 0) / rs_corr[s](i + 1, 0));
                // if(tmp > 1.01)
                if ((tmp > 1. && i < Nt / 2 - 1) || (tmp < 1. && i > Nt / 2 + 1))
                {
                    rs_meff[s](i) = fabs(solver.solve(CoshMeffHelper(
                                                          i - Nt / 2.,
                                                          i + 1 - Nt / 2.,
                                                          tmp
                                                      ), 0, 15).value());
                }
                else
                {
                    rs_meff[s](i) = log(tmp);
                }
            }
        }
    }
    else if (type == MeffType::SINH)
    {
        Roots::BrentRootFinder<double> solver;

        rs_meff.resizeMatrix(Nt - 2, 1);
        for (int s = 0; s < nboot; ++s)
        {
            for (int i = 0; i < Nt - 2; ++i)
            {
                double tmp = fabs(rs_corr[s](i, 0) / rs_corr[s](i + 1, 0));
                if (i < Nt / 2 - 1 && tmp > 1.1 - 1. / (i + 1 - Nt / 2.))
                {
                    rs_meff[s](i) = solver.solve(SinhMeffHelper(
                                                     i - Nt / 2.,
                                                     i + 1 - Nt / 2.,
                                                     tmp
                                                 ), 0.01, 10).value();
                }
                else if (i > Nt / 2 && tmp > 1.1 - 1. / (i + 1 - Nt / 2.))
                {
                    rs_meff[s](i) = solver.solve(SinhMeffHelper(
                                                     i - Nt / 2.,
                                                     i + 1 - Nt / 2.,
                                                     tmp
                                                 ), -10, -0.01).value();
                }
                else
                {
                    rs_meff[s](i) = log(tmp);
                }
            }
        }
    }
    else if (type == MeffType::LOG)
    {
        rs_meff.resizeMatrix(Nt - 2, 1);
        for (int s = 0; s < nboot; ++s)
        {
            for (int i = 0; i < Nt - 2; ++i)
            {
                rs_meff[s](i) = (std::log(fabs(rs_corr[s](i, 0) / rs_corr[s](i + 1, 0))));
            }
        }
    }

    return rs_meff;
}

Sample<Matrix<double>> gevpMeff(
                        MeffType type,
                        const Sample<Matrix<complex<double>>> &rs_gev, int Nt, int t0)
{
    int nboot = rs_gev.size();
    int ngev = rs_gev[0].cols();
    Sample<Matrix<double>> rs_meff(nboot);

    if (type == MeffType::COSH)
    {
        Roots::BrentRootFinder<double> solver;

        rs_meff.resizeMatrix(Nt / 2, ngev);
        for (int s = 0; s < nboot; ++s)
        {
            for (int i = 0; i < Nt / 2; ++i)
            {
                for (int j = 0; j < ngev; ++j)
                {
                    double tmp = fabs(real(rs_gev[s](i, j)));
                    if (tmp > 1.01)
                    {
                        rs_meff[s](i, j) = solver.solve(CoshMeffHelper(
                                                            i - Nt / 2.,
                                                            t0 - Nt / 2.,
                                                            tmp
                                                        ), 0., 10.).value();
                    }
                    else
                    {
                        rs_meff[s](i, j) = log(tmp / fabs(real(rs_gev[s](i + 1, j))));
                    }
                }
            }
        }
    }
    else if (type == MeffType::LOG)
    {
        rs_meff.resizeMatrix(Nt / 2, ngev);
        for (int s = 0; s < nboot; ++s)
        {
            for (int i = 0; i < Nt / 2; ++i)
            {
                for (int j = 0; j < ngev; ++j)
                {
                    rs_meff[s](i, j) = std::log(fabs(real(rs_gev[s](i, j)))) / (t0 - i);
                    // rs_meff[s](i, j) = log(fabs(real(rs_gev[s](i, j)) / real(rs_gev[s](i + 1, j))));
                }
            }
        }
    }

    return rs_meff;
}

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

plat_range findPlateau(const XYDataInterface<double> &data)
{
    const double p_val_threshold = 0.2;

    plat_range res {0, 0};

    index_t start = 0;
    index_t len = 3;

    Models::Line *line = new Models::Line();
    Models::Constant *constant = new Models::Constant();
    vector<double> p_init = {0., 0.};

    Chi2Fit<double, MIN::MIGRAD> Fit(data);
    Fit.options.verbosity = SILENT;

    index_t best_start = start;
    double p_val_min = 1;
    for (; start < data.nPoints() - len; ++start)
    {
        Fit.fitAllPoints(false);
        Fit.fitPointRange(start, start + len - 1, true);
        auto fit = Fit.fit(*constant, p_init);

        double chi2_dof = fit.cost() / fit.nDOF();
        // cout << "chi2_dof = " << chi2_dof << endl;
        double x_mean = mean(data.x({start, start + len}, 0));
        double x_var = variance(data.x({start, start + len}, 0));
        // cout << "x_var = " << x_var << endl;
        double p0_p_se = sqrt(chi2_dof * (1. / len + x_mean * x_mean / x_var));
        // double p1_p_se = fit.err(0);
        double p0_p_val = pvalue(fabs(fit.p(0)), 0., p0_p_se, fit.nDOF());

        cout << "start = " << start << endl;
        cout << "p0_p_se = " << p0_p_se << endl;
        cout << "p_val = " << p0_p_val << endl;

        if (p0_p_val < p_val_min)
        {
            p_val_min = p0_p_val;
            best_start = start;
        }
    }

    res.start = best_start;

    len = data.nPoints() - best_start;
    index_t best_len = len;
    for (; len > 3; --len)
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

        if (p1_p_val > p_val_threshold)
            best_len = len;
    }

    res.len = best_len;

    return res;
}

Sample<double> fitPlateau(
    const LQCDA::Sample<LQCDA::Matrix<double>> &rs_meff,
    const fit_range &f_range)
{
    int nboot = rs_meff.size();
    int npts = rs_meff[0].size();

    Sample<double> result(nboot);

    Vector<double> tvec(npts);
    for (int t = 0; t < npts; ++t)
    {
        tvec(t) = t;
    }

    // fill XYDataSample
    XYDataSample<double> *meff_t = new XYDataSample<double>(npts, 1, 1, nboot);
    for (int n = 0; n < nboot; ++n)
    {
        meff_t->x({}, 0)[n] << tvec;
        meff_t->y({}, 0)[n] << rs_meff[n];
    }
    meff_t->setCovFromSample();

    // search plateau
    plat_range range;
    if (f_range.tmin <= 0 && f_range.tmax <= 0)
    {
        cout << "searching plateau...\n";
        XYData<double> *meff_t_mean = new XYData<double>(npts, 1, 1);
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


    Models::Constant *model = new Models::Constant();

    Chi2Fit<double, MIN::MIGRAD> Fit;
    Fit.options.verbosity = SILENT;

    Sample<double> chi2_dof(nboot);

    vector<double> E0 = {0.};
    for (int n = 0; n < nboot; ++n)
    {
        // Fit E
        Fit.setData(meff_t->getData(n));
        Fit.fitPointRange(range.start, range.start + range.len - 1);

        auto fit = Fit.fit(*model, E0);
        chi2_dof[n] = fit.cost() / fit.nDOF();

        // // Store results
        double E = fabs(fit.parameters()[0]);
        result[n] = E;
    }
    cout << "chi2/dof = " << chi2_dof.mean() << endl;

    // write chi2
    ofstream chi2ofile("chi2.boot");
    FOR_SAMPLE(chi2_dof, s)
    {
        chi2ofile << chi2_dof[s] << endl;
    }

    delete meff_t;
    delete model;

    return result;
}

Sample<double> fitPlateau(
    const LQCDA::Sample<LQCDA::Matrix<double>> &rs_meff,
    const fit_range &f_range,
    Vector<bool> &is_valid)
{
    int nboot = rs_meff.size();
    int npts = rs_meff[0].size();

    Sample<double> result(nboot);

    Vector<double> tvec(npts);
    for (int t = 0; t < npts; ++t)
    {
        tvec(t) = t;
    }

    // fill XYDataSample
    XYDataSample<double> *meff_t = new XYDataSample<double>(npts, 1, 1, nboot);
    for (int n = 0; n < nboot; ++n)
    {
        meff_t->x({}, 0)[n] << tvec;
        meff_t->y({}, 0)[n] << rs_meff[n];
    }
    meff_t->setCovFromSample();

    // search plateau
    plat_range range;
    if (f_range.tmin <= 0 && f_range.tmax <= 0)
    {
        cout << "searching plateau...\n";
        XYData<double> *meff_t_mean = new XYData<double>(npts, 1, 1);
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


    Models::Constant *model = new Models::Constant();

    Chi2Fit<double, MIN::MIGRAD> Fit;
    Fit.options.verbosity = SILENT;

    double chi2_dof = 0.;

    vector<double> E0 = {0.};
    for (int n = 0; n < nboot; ++n)
    {
        // // Fit E
        Fit.setData(meff_t->getData(n));
        Fit.fitPointRange(range.start, range.start + range.len - 1);

        auto fit = Fit.fit(*model, E0);
        chi2_dof += fit.cost() / fit.nDOF();

        // // Store results
        double E = fabs(fit.parameters()[0]);
        double E_fit_err = fabs(fit.errors()[0]);
        result[n] = E;
        is_valid[n] = (E_fit_err / E < 0.5);

    }
    chi2_dof /= nboot;
    cout << "chi2/dof = " << chi2_dof << endl;

    delete meff_t;
    delete model;

    return result;
}

Sample<Matrix<double>> fitCorrelator(
                        const LQCDA::Sample<LQCDA::Matrix<double>> &rs_corr,
                        CorrType type,
                        MeffType mefftype,
                        unsigned int nStates,
                        const fit_range &f_range)
{
    int nboot = rs_corr.size();
    int npts = rs_corr[0].rows();

    Sample<Matrix<double>> result(nboot, 2, nStates);

    Vector<double> tvec(npts);
    for (int t = 0; t < npts; ++t)
    {
        tvec(t) = t;
    }

    // fill XYDataSample
    XYDataSample<double> *corr_t = new XYDataSample<double>(npts, 1, 1, nboot);
    for (int n = 0; n < nboot; ++n)
    {
        corr_t->x({}, 0)[n] << tvec;
        corr_t->y({}, 0)[n] << rs_corr[n].col(0);
    }
    corr_t->setCovFromSample();

    cout << "using plateau range (" << f_range.tmin << ", "
             << f_range.tmax << ")\n";

    // fit correlator
    if (type == CorrType::EXP)
    {
        Models::multiAExpmEt *model = new Models::multiAExpmEt(nStates);
        Chi2Fit<double, MIN::MIGRAD> Fit;
        Fit.options.verbosity = SILENT;

        double chi2_dof = 0.;

        // get estimates of the energies from meff fit
        auto rs_meff = localMeff(mefftype, rs_corr, npts);
        auto mass = fitPlateau(rs_meff, f_range).mean();

        vector<double> a_E_init(2 * nStates);
        for(int j=0; j<nStates; ++j)
        {
            a_E_init[2*j] = rs_corr.mean()(npts/4, 0)*exp(mass*(npts/4));
            a_E_init[2*j+1] = mass;
        }

        for (int n = 0; n < nboot; ++n)
        {
            // Fit a's and E's
            Fit.setData(corr_t->getData(n));
            Fit.fitPointRange(f_range.tmin, f_range.tmax);

            auto fit = Fit.fit(*model, a_E_init);
            chi2_dof += fit.cost() / fit.nDOF();

            // Store results
            for (unsigned int i = 0; i < nStates; ++i)
            {
                double ai = fabs(fit.parameters()[2*i]);
                double ai_fit_err = fabs(fit.errors()[2*i]);
                double Ei = fabs(fit.parameters()[2*i+1]);
                double Ei_fit_err = fabs(fit.errors()[2*i+1]);
                result[n](0, i) = ai;
                result[n](1, i) = Ei;
            }
        }
        chi2_dof /= nboot;
        cout << "chi2/dof = " << chi2_dof << endl;
        delete model;
    }
    else if (type == CorrType::VARPRO)
    {
        double chi2_dof = 0.;

        // get estimates of the energies from meff fit
        auto rs_meff = localMeff(mefftype, rs_corr, npts);
        auto rs_mass = fitPlateau(rs_meff, automatic_plat_range(3.8, 40, false, false));

        // create minimizer
        auto Minimizer = std::unique_ptr<MIN::DIFFERENTIAL_EVOLUTION<double>>(new MIN::DIFFERENTIAL_EVOLUTION<double>);
        Minimizer->options().verbosity = SILENT;
        Minimizer->dither = false;
        Minimizer->tolerance = 1.e-8;
        Minimizer->max_iter = 1e4;

        vector<double> E_init(nStates);
        vector<double> E_err_init(nStates);
        E_init[0] = rs_mass.mean();
        E_err_init[0] = sqrt(rs_mass.variance());
        for(int j=1; j<nStates; ++j)
        {
            E_init[j] = (1+j*rs_mass.mean())*rs_mass.mean();
            E_err_init[j] = j;
        }

        for (int n = 0; n < nboot; ++n)
        {
            unsigned int nFitPts = f_range.tmax - f_range.tmin + 1;

            // create varpro function
            VarProChi2 *chi2 = new VarProChi2(nStates, f_range, corr_t->getData(n), corr_t->yyCov(0, 0).block(f_range.tmin, f_range.tmin, nFitPts, nFitPts));

            // Fit E's
            auto min = Minimizer->minimize(*chi2, E_init, E_err_init);
            // Store results
            Vector<double> E(nStates);
            for (unsigned int i = 0; i < nStates; ++i)
            {
                E(i) = fabs(min.minimum[i]);
            }

            // Get a's
            Matrix<double> phi(nFitPts, nStates);
            FOR_MAT(phi, i, j)
            {
                phi(i, j) = exp(-E(j)*((*corr_t)[n].x(i, 0)));
            }
            Vector<double> a = PseudoInverse(phi) * (*corr_t)[n].y({f_range.tmin, f_range.tmax+1}, 0);

            // Store results
            for (unsigned int i = 0; i < nStates; ++i)
            {
                result[n](0, i) = a(i);
                result[n](1, i) = E(i);
            }

            double nDOF = corr_t->yDim() * nFitPts - nStates;
            chi2_dof += min.final_cost / nDOF;
        }
        chi2_dof /= nboot;
        cout << "chi2/dof = " << chi2_dof << endl;
    }

    delete corr_t;

    return result;
}

Sample<double> getPlateau(
    const Sample<Matrix<double>> &rs_corr,
    MeffType type,
    const fit_range &f_range)
{
    int nboot = rs_corr.size();
    int T = rs_corr[0].rows();

    Sample<Matrix<double>> rs_meff = localMeff(type, rs_corr, T);

    int npts = rs_meff[0].size();

    Sample<double> result(nboot);

    Vector<double> tvec(npts);
    for (int t = 0; t < npts; ++t)
    {
        tvec(t) = t;
    }

    // fill XYDataSample
    XYDataSample<double> *meff_t = new XYDataSample<double>(npts, 1, 1, nboot);
    for (int n = 0; n < nboot; ++n)
    {
        meff_t->x({}, 0)[n] << tvec;
        meff_t->y({}, 0)[n] << rs_meff[n];
    }
    meff_t->setCovFromSample();

    // search plateau
    plat_range range;
    if (f_range.tmin <= 0 && f_range.tmax <= 0)
    {
        cout << "searching plateau...\n";
        XYData<double> *meff_t_mean = new XYData<double>(npts, 1, 1);
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


    Models::Constant *model = new Models::Constant();

    Chi2Fit<double, MIN::MIGRAD> Fit;
    Fit.options.verbosity = SILENT;

    double chi2_dof = 0.;

    vector<double> E0 = {0.};
    for (int n = 0; n < nboot; ++n)
    {
        // // Fit E
        Fit.setData(meff_t->getData(n));
        Fit.fitPointRange(range.start, range.start + range.len - 1);

        auto fit = Fit.fit(*model, E0);
        chi2_dof += fit.cost() / fit.nDOF();

        // // Store results
        double E = fabs(fit.parameters()[0]);
        double E_fit_err = fabs(fit.errors()[0]);
        result[n] = E;

    }
    chi2_dof /= nboot;
    cout << "chi2/dof = " << chi2_dof << endl;

    delete meff_t;
    delete model;

    return result;
}

#endif // PLATEAU_CPP