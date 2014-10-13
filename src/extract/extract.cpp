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

using namespace std;
using namespace LQCDA;
namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;
namespace fs = boost::filesystem;

bool has_t_column(const Sample<Matrix<double>>& sample)
{
    bool res = true;
    Matrix<double> first_col = sample.mean().col(0);
    for(int i=0; i<sample.rows(); i++)
    {
        res &= (first_col(i,0) == i);
    }
    return res;
}

int extract(const std::string& manfile, const extract_parameters& params)
{
    // Print info
    cout << endl << endl
    << "********************************************************************" << endl
    << "Extracting:" << endl
    << "Manifest file " << manfile << endl
    << "Interpolator " << params.corr_header << endl
    << "********************************************************************" << endl;

    fs::path manfilepath(manfile);

    // Get files from manfile
    vector<string> prop_files;
    ifstream man_if(manfile);
    string buf;
    while(getline(man_if, buf))
    {
        prop_files.push_back(buf);
    }

    // Read correlators for all configurations
    Sample<Matrix<double>> corr_sample;
    readSamples(&corr_sample, manfile, params.corr_header);

    if(has_t_column(corr_sample))
        corr_sample = corr_sample.block(0, 1, corr_sample.rows(), corr_sample.cols()-1);

    // Generate random number generator state
    RandGen rng;
    RandGen::rg_state state;
    rng.getState(state);

    // Resample
    Sample<Matrix<double>> rs_corr = resample(corr_sample, params.nboot, state);

    // Print correlator
    cout << endl << endl
    << "********************************************************************" << endl
    << "Correlator" << endl;

    ostream * corr_out;
    if(params.save_corr)
    {
        ofstream *  fout = new ofstream(params.corr_file);
        tee_device<ostream, ofstream> * tee_out = new tee_device<ostream, ofstream>(cout, *fout);
        corr_out = new stream<tee_device<ostream, ofstream>>(*tee_out);
    }
    else
    {
        corr_out = new ostream(cout.rdbuf());
    }

    Matrix<double> corr = rs_corr.mean();
    Matrix<double> corr_var = rs_corr.variance();
    for(int i=0; i<corr.rows(); i++)
    {
        (*corr_out) << i;
        for(int j=0; j<corr.cols(); j++)
        {
            (*corr_out) << ' ' << corr(i, j) << ' ' << sqrt(corr_var(i, j));
        }
        (*corr_out) << endl;
    }

    unsigned int Nt = corr.rows();

    // Extract effective mass if asked
    if(params.extract_meff)
    {
        // Compute effective mass
        Sample<Matrix<double>> rs_meff = localMeff(params.meff_type, rs_corr, Nt);

        // Print effective mass
        cout << endl << endl
        << "********************************************************************" << endl
        << "Effective mass" << endl;
        ostream * meff_out;
        if(params.save_meff)
        {
            ofstream *  fout = new ofstream(params.meff_file);
            tee_device<ostream, ofstream> * tee_out = new tee_device<ostream, ofstream>(cout, *fout);
            meff_out = new stream<tee_device<ostream, ofstream>>(*tee_out);
        }
        else
        {
            meff_out = new ostream(cout.rdbuf());
        }

        Matrix<double> meff = rs_meff.mean();
        Matrix<double> meff_var = rs_meff.variance();
        for(int i=0; i<meff.rows(); i++)
        {
            (*meff_out) << i;
            for(int j=0; j<meff.cols(); j++)
            {
                (*meff_out) << ' ' << meff(i, j) << ' ' << sqrt(meff_var(i, j));
            }
            (*meff_out) << endl;
        }
    }

    return 0;
}