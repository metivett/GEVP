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
	// number of bootstrap samples
 	unsigned int nboot;
 	// pion physical mass (GeV)
 	double pi_mass;
 	// K physical mass (GeV)
 	double K_mass;
 	// Omega physical mass (GeV)
 	double Omega_mass;
};

std::ostream& operator<<(std::ostream& os, const set_scale_parameters& p);

struct set_scale_result
{
	double a;
	double a_err;
};

namespace Models {
    class aMOmega_vs_aMpi_aMKchi
    : public LQCDA::ParametrizedScalarFunction<double>
    {
    public:
        aMOmega_vs_aMpi_aMKchi(double Mpi_phys, double MK_phys, double MOmega_phys)
        : LQCDA::ParametrizedScalarFunction<double>(2, 3)
        , m_Mpi_phys(Mpi_phys)
        , m_MKchi2_phys(MK_phys*MK_phys - Mpi_phys*Mpi_phys/2.)
        , m_MOmega_phys(MOmega_phys)
        {
            double MOmega2_phys =  m_MOmega_phys*m_MOmega_phys;
            m_Mpi2_MOmega2_phys = m_Mpi_phys*m_Mpi_phys / MOmega2_phys;
            m_MKchi2_MOmega2_phys = m_MKchi2_phys / MOmega2_phys;

        }

        virtual double operator()(const double* x, const double* p) const override
        {
            return  p[0]*m_MOmega_phys 
                    + p[1]*(x[0] - m_Mpi2_MOmega2_phys)
                    + p[2]*(x[1] - m_MKchi2_MOmega2_phys)
            ;
        }

    private:
        double m_Mpi_phys;
        double m_MKchi2_phys;
        double m_MOmega_phys;
        double m_Mpi2_MOmega2_phys;
        double m_MKchi2_MOmega2_phys;
    };
}

set_scale_result set_scale(const std::string& manifest, const set_scale_parameters& parameters);

#endif // SET_SCALE_HPP
