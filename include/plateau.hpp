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

namespace Models
{

class Constant
    : public LQCDA::ParametrizedScalarFunction<double>
{
public:
    Constant()
        : LQCDA::ParametrizedScalarFunction<double>(1, 1)
    {}

    virtual double operator()(const double *x, const double *p) const override
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

    virtual double operator()(const double *x, const double *p) const override
    {
        return p[0] + p[1] * (*x);
    }
};

class multiAExpmEt
    : public LQCDA::ParametrizedScalarFunction<double>
{
public:
    multiAExpmEt(unsigned int nExp)
        : LQCDA::ParametrizedScalarFunction<double>(1, 2 * nExp)
        , m_nExp(nExp)
    {}

    virtual double operator()(const double *x, const double *p) const override
    {
        double res = 0.;
        for (unsigned int n = 0; n < m_nExp; ++n)
        {
            res += p[2 * n] * exp(-p[2 * n + 1] * x[0]);
        }
        return res;
    }
private:
    unsigned m_nExp;
};

}

class VarProChi2
    : public LQCDA::ScalarFunction<double>
{
public:
    VarProChi2(unsigned int nExp, const fit_range &f_range, const LQCDA::XYDataInterface<double> &data, LQCDA::ConstRef<LQCDA::Matrix<double>> cov)
        : LQCDA::ScalarFunction<double>(nExp)
        , m_nExp(nExp)
        , m_fRange(f_range)
        , m_nPts(f_range.tmax - f_range.tmin + 1)
        , m_Data(&data)
        , m_Cinv(cov.inverse())
    {}

    virtual double operator()(const double *E) const override
    {
        LQCDA::Matrix<double> phi(m_nPts, m_nExp);
        for (unsigned int j = 0; j < phi.cols(); ++j)
            for (unsigned int i = 0; i < phi.rows(); ++i)
            {
                phi(i, j) = exp(-E[j]*m_Data->x(i, 0));
            }

        LQCDA::Vector<double> r = m_Data->y({m_fRange.tmin, m_fRange.tmax+1}, 0) - phi*LQCDA::PseudoInverse(phi)*m_Data->y({m_fRange.tmin, m_fRange.tmax+1}, 0);

        double chi2 = r.dot(m_Cinv * r);
        return chi2;
    }

private:
    unsigned int m_nExp;
    fit_range m_fRange;
    unsigned int m_nPts;
    const LQCDA::XYDataInterface<double> *m_Data;
    // inverse covariance matrix
    LQCDA::Matrix<double> m_Cinv;
};

enum class MeffType
{
    LOG,
    COSH,
    SINH
};

enum class CorrType
{
    EXP,
    VARPRO
};

fit_range automatic_plat_range(double beta, unsigned int T, bool is_strange, bool delayed = false);

// local effective mass
LQCDA::Sample<LQCDA::Matrix<double>> localMeff(
                                      MeffType type,
                                      const LQCDA::Sample<LQCDA::Matrix<double>> &rs_corr, int Nt);
//gevp effective mass
LQCDA::Sample<LQCDA::Matrix<double>> gevpMeff(
                                      MeffType type,
                                      const LQCDA::Sample<LQCDA::Matrix<std::complex<double>>> &rs_gev, int Nt, int t0);

// fit plateau
LQCDA::Sample<double> fitPlateau(
    const LQCDA::Sample<LQCDA::Matrix<double>> &rs_meff,
    const fit_range &range = fit_range());

LQCDA::Sample<double> fitPlateau(
    const LQCDA::Sample<LQCDA::Matrix<double>> &rs_meff,
    const fit_range &f_range,
    LQCDA::Vector<bool> &is_valid);

// fit correlator
LQCDA::Sample<LQCDA::Matrix<double>> fitCorrelator(
                                      const LQCDA::Sample<LQCDA::Matrix<double>> &rs_corr,
                                      CorrType type,
                                      MeffType mefftype,
                                      unsigned int nStates,
                                      const fit_range &f_range);

// get plateau (from correlator)
LQCDA::Sample<double> getPlateau(
    const LQCDA::Sample<LQCDA::Matrix<double>> &rs_corr,
    MeffType type,
    const fit_range &range = fit_range());


#endif // PLATEAU_HPP