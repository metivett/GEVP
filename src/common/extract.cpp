/*
 * extract.cpp
 *
 *  Created on: Sep 15, 2014
 *      Author: Thibaut Metivet
 */

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/filesystem.hpp>

#include "extract.hpp"
#include "hadron.hpp"

using namespace std;
using namespace LQCDA;
namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;
namespace fs = boost::filesystem;

bool has_t_column(const Sample<Matrix<double>> &sample)
{
    bool res = true;
    Matrix<double> first_col = sample.mean().col(0);
    for (int i = 0; i < sample.rows(); i++)
    {
        res &= (first_col(i, 0) == i);
    }
    return res;
}

vector<string> get_files_from_manfile(const std::string &manfile)
{
    vector<string> prop_files;
    ifstream man_if(manfile);
    string buf;
    while (getline(man_if, buf))
    {
        prop_files.push_back(buf);
    }
    return prop_files;
}

std::string header_from_particle(Hadron part, const std::string &smearing)
{
    using Part = Hadron;
    std::string op;
    switch (part)
    {
    case Part::pion: op = "PP00"; break;
    case Part::kaon: op = "PP01"; break;
    case Part::rho: op = "V1V100"; break;
    case Part::nucleon: op = "NUCL00"; break;
    case Part::sigma: op = "NUCL10"; break;
    case Part::xi: op = "NUCL01"; break;
    case Part::omega: op = "DELTA11"; break;
    case Part::xi_star: op = "DELTA01"; break;
    case Part::sigma_star: op = "DELTA10"; break;
    case Part::delta: op = "DELTA00"; break;
    }
    op += smearing;
    return op;
}

Sample<Matrix<double>> read_corr_sample(const std::string &manfile, Hadron part, const std::string &smearing)
{
    Sample<Matrix<double>> corr_sample;
    readSamples(&corr_sample, manfile, header_from_particle(part, smearing));

    // if (hadron_properties[part].is_baryon)
    if(false) // Do NOT project on parity +1
    {
        Sample<Matrix<double>> parity_corr_sample(corr_sample.size(), corr_sample.rows(), std::ceil(corr_sample.cols() / 2.));
        FOR_SAMPLE(corr_sample, s)
        {
            if (has_t_column(corr_sample))
            {
                parity_corr_sample[s].col(0) = corr_sample[s].col(0);
                parity_corr_sample[s].col(1) = corr_sample[s].col(1) - corr_sample[s].col(3);
                parity_corr_sample[s].col(2) = corr_sample[s].col(2) - corr_sample[s].col(4);
            }
            else
            {
                parity_corr_sample[s].col(0) = corr_sample[s].col(0) - corr_sample[s].col(2);
                parity_corr_sample[s].col(1) = corr_sample[s].col(1) - corr_sample[s].col(3);
            }
        }
        return parity_corr_sample;
    }
    else
    {
        return corr_sample;
    }
}

Hadron hadron_from_string(const string &s)
{
    using Part = Hadron;
    if (s == "pion") return Part::pion;
    else if (s == "kaon") return Part::kaon;
    else if (s == "rho") return Part::rho;
    else if (s == "nucleon") return Part::nucleon;
    else if (s == "sigma") return Part::sigma;
    else if (s == "xi") return Part::xi;
    else if (s == "omega") return Part::omega;
    else if (s == "xi_star") return Part::xi_star;
    else if (s == "sigma_star") return Part::sigma_star;
    else if (s == "delta") return Part::delta;
    else throw std::runtime_error("bad particle type");
}

Sample<Matrix<double>> extract_corr(const std::string &manfile, const extract_parameters &params)
{
    fs::path manfilepath(manfile);

    // Get files from manfile
    vector<string> prop_files = get_files_from_manfile(manfile);

    // Read correlators for all configurations
    Sample<Matrix<double>> corr_sample = read_corr_sample(manfile, params.particle, params.smearing);

    if (has_t_column(corr_sample))
        corr_sample = corr_sample.block(0, 1, corr_sample.rows(), corr_sample.cols() - 1);

    // Generate random number generator state
    RandGen rng;
    RandGen::rg_state state;
    rng.getState(state);

    // Resample
    Sample<Matrix<double>> rs_corr = resample(corr_sample, params.nboot, state);

    unsigned int Nt = rs_corr[0].rows();

    // Fold if requested
    if(params.fold_correlator)
    {
        FOR_SAMPLE(rs_corr, s)
        {
            for(int t = 1; t < Nt/2; ++t)
            {
                rs_corr[s].row(t) = (rs_corr[s].row(t) + rs_corr[s].row(Nt-t))/2.;

            }
        }
    }

    return rs_corr;
}

Sample<Matrix<double>> extract_meff(const std::string &manfile, const extract_parameters &params, MeffType mefftype)
{
    Sample<Matrix<double>> rs_corr = extract_corr(manfile, params);
    unsigned int Nt = rs_corr[0].rows();

    // Compute effective mass
    Sample<Matrix<double>> rs_meff = localMeff(mefftype, rs_corr, Nt);

    return rs_meff;
}

Sample<double> extract_mass(const std::string &manfile, const extract_parameters &params, MeffType mefftype, const fit_range &range)
{
    Sample<Matrix<double>> rs_meff = extract_meff(manfile, params, mefftype);

    fit_range f_range;
    if(params.beta > 0. && (range.tmin <= 0 && range.tmax <= 0))
        f_range = automatic_plat_range(params.beta, rs_meff.rows(), hadron_properties[params.particle].is_strange, params.delay);
    else 
        f_range = range;

    Sample<double> rs_m = fitPlateau(rs_meff, f_range);
    return rs_m;
}

Sample<Matrix<double>> extract_mass_varpro(const std::string &manfile, const extract_parameters &params, MeffType mefftype, unsigned int nExp, const fit_range &range)
{
    Sample<Matrix<double>> rs_corr = extract_corr(manfile, params);

    fit_range f_range;
    if(params.beta > 0. && (range.tmin <= 0 && range.tmax <= 0))
        f_range = automatic_plat_range(params.beta, rs_corr.rows(), hadron_properties[params.particle].is_strange, params.delay);
    else 
        f_range = range;

    Sample<Matrix<double>> rs_a_E = fitCorrelator(rs_corr, CorrType::VARPRO, mefftype, nExp, f_range);

    return rs_a_E;
}