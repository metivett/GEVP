/*
 * gevp.cpp
 *
 *  Created on: Apr 24, 2014
 *      Author: Thibaut Metivet
 */

#include "gevp.hpp"
#include "Z001.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <complex>

using namespace LQCDA;

// put a line of reals into the gevp form (square complex matrix)
template<typename Derived>
Matrix<std::complex<double>> gevp_form(const MatrixExpr<Derived> &mat)
{
    int rows = sqrt(mat.cols() / 2);
    Matrix<std::complex<double>> res(rows, rows);
    FOR_MAT(res, i, j)
    {
        index_t n = i * rows + j;
        res(i, j) = std::complex<double>(mat(0, 2 * n), mat(0, 2 * n + 1));
    }

    return res;
}

double real(const std::complex<double> &c)
{
    return c.real();
}

Matrix<std::complex<double>> gev(const Matrix<double> &mat, unsigned int t0)
{
    int Nt = mat.rows();
    unsigned int ngev = sqrt(mat.cols() / 2);

    Matrix<std::complex<double>> result(Nt, ngev);

    Matrix<std::complex<double>> Ct0_inv = gevp_form(mat.row(t0)).inverse();
    Matrix<std::complex<double>> Ct;

    for (int t = 0; t < Nt; ++t)
    {
        Ct = gevp_form(mat.row(t));
        Eigen::ComplexEigenSolver<Matrix<std::complex<double>>> eigsolv(Ct0_inv * Ct, false);
        result.row(t) = eigsolv.eigenvalues().transpose();
        // std::cout << Ct << '\n';
    }

    return result;
}

Matrix<double> gev_svd(const Matrix<double> &mat, unsigned int t0)
{
    int Nt = mat.rows();
    unsigned int ngev = sqrt(mat.cols() / 2);

    Matrix<double> result(Nt, ngev);

    Matrix<std::complex<double>> Ct0_inv = gevp_form(mat.row(t0)).inverse();
    Matrix<std::complex<double>> Ct;

    for (int t = 0; t < Nt; ++t)
    {
        Ct = gevp_form(mat.row(t));
        // Eigen::ComplexEigenSolver<Matrix<std::complex<double>>> eigsolv(Ct0_inv * Ct, false);
        Eigen::JacobiSVD<Matrix<std::complex<double>>> svd(Ct0_inv * Ct);
        // result.row(t) = eigsolv.eigenvalues().transpose();
        result.row(t) = svd.singularValues().transpose();
        // std::cout << Ct << '\n';
    }

    return result;
}

double deltaM2_COM(double W, double MPi, double L)
{
    double E_COM = W;
    double q2 = (E_COM * E_COM / 4. - MPi * MPi) * (L * L) / (4.*M_PI * M_PI);
    double delta;
    if (q2 < 0.)
    {
        delta = 0.;
        std::cout << "q2 < 0 !" << std::endl;
    }
    else
        delta = -8.*M_PI * M_SQRTPI / (L * L * L) * q2 / E_COM * Z001::z001q2(q2);
    // std::cout << "q2 = " << q2 << std::endl;
    return delta;
}

// Luscher's cotg(delta)
double cotd_COM(double W, double MPi, double L)
{
    double E_COM = W;
    double q2 = (E_COM * E_COM / 4. - MPi * MPi) * (L * L) / (4.*M_PI * M_PI);
    double cotd;
    if (q2 < 0.)
    {
        cotd = 0.;
        std::cout << "q2 < 0 !" << std::endl;
    }
    else
        cotd = Z001::z001q2(q2) / (M_PI * M_SQRTPI * sqrt(q2));
    // std::cout << "q2 = " << q2 << std::endl;
    return cotd;
}

Sample<Matrix<double>> fitSin2d(
                        const Matrix<double> &gevp_plat,
                        const Sample<Matrix<double>> &rs_sin2d,
                        const Sample<double> &rs_pi_plat,
                        const fit_range &range)
{
    // std::cout << rs_sin2d[0].rows() << std::endl << rs_sin2d[0](4, 0);
    int nboot = rs_sin2d.size();
    int npts = gevp_plat.size();

    Sample<Matrix<double>> result(nboot, 2, 2);

    XYDataSample<double> *sin2d_vs_E_mpi = new XYDataSample<double>(npts, 2, 1, nboot);

    for (int n = 0; n < nboot; ++n)
    {
        sin2d_vs_E_mpi->x({}, 0)[n] << gevp_plat.transpose();
        sin2d_vs_E_mpi->x({}, 1)[n] << Matrix<double>::Constant(npts, 1, rs_pi_plat[n]);
        sin2d_vs_E_mpi->y({}, 0)[n] << rs_sin2d[n];
    }
    // sin2d_vs_E_mpi->setCovFromSample();
    Matrix<double> yycov(npts, npts);
    yycov.setZero();
    auto sin2dvar = rs_sin2d.variance();
    FOR_MAT_DIAG(yycov, i)
    {
        yycov(i, i) = sin2dvar(i, 0);
    }
    sin2d_vs_E_mpi->yyCov(0, 0) = yycov;
    std::cout << yycov << std::endl;

    Models::Sin2d_vs_E_mpi *model = new Models::Sin2d_vs_E_mpi();

    Chi2Fit<double, MIN::MIGRAD> Fit;
    Fit.options.verbosity = SILENT;

    double chi2_dof = 0.;

    std::vector<double> g2_mrho2_init = {36, 0.22};
    // std::vector<double> g2_mrho2_init = {20., 0.23};
    // std::vector<double> g2_mrho2_init = {36., 0.075};
    // std::vector<double> g2_mrho2_init = {36., 0.275};
    // std::vector<double> g2_mrho2_init = {60., 0.119};

    // Test: fit mean
    Chi2Fit<double, MIN::MIGRAD> MeanFit;
    MeanFit.options.verbosity = SILENT;
    auto meandata = sin2d_vs_E_mpi->mean();
    MeanFit.setData(meandata);
    MeanFit.fitPointRange(range.tmin, range.tmax, true);
    auto meanfit = MeanFit.fit(*model, g2_mrho2_init);
    std::cout << "mean fit: " << std::endl
              << "g2=" << meanfit.parameters()[0] 
              << " mrho2=" << meanfit.parameters()[1]<< std::endl;

    std::vector<double> g2_mrho2_from_mean = {fabs(meanfit.parameters()[0]), fabs(meanfit.parameters()[1])};

    for (int n = 0; n < nboot; ++n)
    {
        // Fit g and Mrho
        Fit.setData(sin2d_vs_E_mpi->getData(n));

        if (range.tmin <= 0 && range.tmax <= 0)
            Fit.fitAllPoints(true);
        else
            Fit.fitPointRange(range.tmin, range.tmax, true);

        auto fit = Fit.fit(*model, g2_mrho2_from_mean);

        // Store results
        if (fit.isValid())
        {
            double g2 = fabs(fit.parameters()[0]);
            double Mrho2 = fabs(fit.parameters()[1]);
            result[n](0, 0) = g2;
            result[n](0, 1) = Mrho2;
            result[n](1, 0) = fabs(fit.errors()[0]);
            result[n](1, 1) = fabs(fit.errors()[1]);

            chi2_dof += fit.cost() / fit.nDOF();
        }
        else
        {
            ScalarConstraint<double> fixedValue;
            // Fix Mrho2 with g2 fixed
            fixedValue.fixValue(g2_mrho2_from_mean[0]);
            auto fitMrho2 = Fit.fit(*model, {36, fabs(fit.parameters()[1])}, {fixedValue, ScalarConstraint<double>()});

            std::cout << fitMrho2.isValid() << " mrho2=" << fitMrho2.parameters()[1] << " ";

            // Fit g2 with obtained value for Mrho2 fixed
            fixedValue.fixValue(fabs(fitMrho2.parameters()[1]));
            auto fitg2 = Fit.fit(*model, {fabs(fit.parameters()[0]), fabs(g2_mrho2_from_mean[1])}, {ScalarConstraint<double>(), fixedValue});

            std::cout << fitg2.isValid() << " g2=" << fitg2.parameters()[0] << " ";

            // Complete fit with obtained initial values
            std::vector<double> g2_mrho2_init2 = {fabs(fitg2.parameters()[0]), fabs(fitMrho2.parameters()[1])};
            auto fit2 = Fit.fit(*model, g2_mrho2_init2);

            std::cout << fit2.isValid() << std::endl;
            // std::cout << fit2.isValid() << "\t" << sqrt(fabs(fit2.parameters()[0])) << "\t" << fabs(fit2.errors()[0])/fabs(fit2.parameters()[0]) << std::endl;

            chi2_dof += fit2.cost() / fit2.nDOF();
            if (fit2.isValid())
            {
                result[n](0, 0) = fabs(fit2.parameters()[0]);
                result[n](0, 1) = fabs(fit2.parameters()[1]);
                result[n](1, 0) = fabs(fit2.errors()[0]);
                result[n](1, 1) = fabs(fit2.errors()[1]);
            }
            else
            {
                result[n](0, 0) = fabs(fit2.parameters()[0]);
                result[n](0, 1) = fabs(fit2.parameters()[1]);
                result[n](1, 0) = fabs(fit2.parameters()[0]);
                result[n](1, 1) = fabs(fit2.parameters()[1]);
            }
        }
    }
    chi2_dof /= nboot;
    std::cout << "chi2/dof = " << chi2_dof << std::endl;

    delete sin2d_vs_E_mpi;
    delete model;

    return result;
}