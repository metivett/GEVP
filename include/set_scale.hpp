/*
 * set_scale.hpp
 *
 *  Created on: Jul 21, 2014
 *      Author: Thibaut Metivet
 */

#ifndef SET_SCALE_HPP
#define SET_SCALE_HPP

#include "LQCDA.hpp"
#include "utils.hpp"

struct set_scale_parameters
{
    // pion meson plateau fit range
    fit_range pi_fit_range;
    // K meson plateau fit range
    fit_range K_fit_range;
    // Omega plateau fit range
    fit_range Omega_fit_range;
    // beta
    double beta;
    // number of bootstrap samples
    unsigned int nboot;
    // pion physical mass (GeV)
    double pi_mass;
    // K physical mass (GeV)
    double K_mass;
    // Omega physical mass (GeV)
    double Omega_mass;
};

struct set_multi_scales_parameters
{
    // betas
    std::vector<double> betas;
    // number of bootstrap samples
    unsigned int nboot;
    // pion physical mass (GeV)
    double pi_mass;
    // K physical mass (GeV)
    double K_mass;
    // Omega physical mass (GeV)
    double Omega_mass;
};

std::ostream &operator<<(std::ostream &os, const set_scale_parameters &p);


namespace Models
{
class aMOmega_vs_aMpi2_aMKchi2
    : public LQCDA::ParametrizedScalarFunction<double>
{
public:
    aMOmega_vs_aMpi2_aMKchi2(double Mpi_phys, double MK_phys, double MOmega_phys)
        : LQCDA::ParametrizedScalarFunction<double>(2, 3)
        , m_Mpi_phys(Mpi_phys)
        , m_MKchi2_phys(MK_phys * MK_phys - Mpi_phys *Mpi_phys / 2.)
        , m_MOmega_phys(MOmega_phys)
    {
        double MOmega2_phys =  m_MOmega_phys * m_MOmega_phys;
        m_Mpi2_MOmega2_phys = m_Mpi_phys * m_Mpi_phys / MOmega2_phys;
        m_MKchi2_MOmega2_phys = m_MKchi2_phys / MOmega2_phys;

    }

    virtual double operator()(const double *x, const double *p) const override
    {
        // return  p[0]*m_MOmega_phys *
        //         (1 + p[1]*(x[0] - m_Mpi2_MOmega2_phys) + p[2]*(x[1] - m_MKchi2_MOmega2_phys))
        //         ;
        return p[0] * m_MOmega_phys
               + p[1] * (x[0] - m_Mpi2_MOmega2_phys) 
               + p[2] * (x[1] - m_MKchi2_MOmega2_phys)
        ;
    }

private:
    double m_Mpi_phys;
    double m_MKchi2_phys;
    double m_MOmega_phys;
    double m_Mpi2_MOmega2_phys;
    double m_MKchi2_MOmega2_phys;
};

// Set several scales simultaneously
// with aMOmega = A(beta)*MOmega_phys*(1 + B*(m_Mpi2_MOmega2-m_Mpi2_MOmega2_phys) + C*(MKchi2_MOmega2-MKchi2_MOmega2_phys))
// x in R^3 : {Mpi2_MOmega2, MKchi2_MOmega2, n_beta}
// p in R^(nscales+2): {A(beta_1), ..., A(beta_nscales), B, C}
class aMOmega_vs_aMpi2_aMKchi2_multi
    : public LQCDA::ParametrizedScalarFunction<double>
{
public:
    aMOmega_vs_aMpi2_aMKchi2_multi(double Mpi_phys, double MK_phys, double MOmega_phys, unsigned int nbetas)
        : LQCDA::ParametrizedScalarFunction<double>(3, nbetas + 2)
        , m_Mpi_phys(Mpi_phys)
        , m_MKchi2_phys(MK_phys * MK_phys - Mpi_phys *Mpi_phys / 2.)
        , m_MOmega_phys(MOmega_phys)
        , m_nBetas(nbetas)
    {
        double MOmega2_phys =  m_MOmega_phys * m_MOmega_phys;
        m_Mpi2_MOmega2_phys = m_Mpi_phys * m_Mpi_phys / MOmega2_phys;
        m_MKchi2_MOmega2_phys = m_MKchi2_phys / MOmega2_phys;

    }

    virtual double operator()(const double *x, const double *p) const override
    {
        double A;
        for(unsigned int nbeta=0; nbeta<m_nBetas; nbeta++)
        {
            if(x[2] == static_cast<double>(nbeta))
            {
                A = p[nbeta];
                break;
            }
        }
        double B = p[m_nBetas];
        double C = p[m_nBetas+1];
        return  A*m_MOmega_phys *
                (1 + B*(x[0] - m_Mpi2_MOmega2_phys) + C*(x[1] - m_MKchi2_MOmega2_phys))
                ;
    }

private:
    double m_Mpi_phys;
    double m_MKchi2_phys;
    double m_MOmega_phys;
    double m_Mpi2_MOmega2_phys;
    double m_MKchi2_MOmega2_phys;
    unsigned int m_nBetas;
};
}

LQCDA::Sample<double> set_scale(const std::string &manifest, const set_scale_parameters &parameters);
LQCDA::Sample<LQCDA::Matrix<double>> set_multi_scales(const std::vector<std::string> &manifests, const set_multi_scales_parameters &parameters);

#endif // SET_SCALE_HPP
