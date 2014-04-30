/*
 * GEVP.cpp
 *
 *  Created on: Jun 14, 2013
 *      Author: Thibaut Metivet
 */

#include "LQCDA.hpp"
#include "GEVP.hpp"
#include "GEVPFitModels.hpp"
#include "Z001.hpp"
#include <functional>
#include <cmath>
#include <iostream>

using namespace LQCDA;

template<class T>
struct Evals
{
public:
    T L1, L2;

    Evals();
};

std::vector<std::vector<Eigen::Matrix2cd> > ResampleGEVPCorrelators(const std::vector<std::vector<Eigen::Matrix2cd> >& sample, unsigned int Nboot, RandGen::rg_state state)
{
    RandGen rng;
    rng.setState(state);
    
    std::vector<std::vector<Eigen::Matrix2cd> > result(Nboot);

    unsigned int NSample = sample.size();
    double NSample_inv = 1. / (double)NSample;
    unsigned int index;
    
    for(int n = 0; n < Nboot; ++n) {
	index = rng.getUniformInt(NSample);
	std::vector<Eigen::Matrix2cd> buf(sample[index]);
	for(int i = 1; i < NSample; ++i) {
	    index = rng.getUniformInt(NSample);
	    buf = buf + sample[index];
	}
	buf = buf * NSample_inv;

	result[n] = buf;
    }

    return result;
}

std::vector<GEVPPlateau> ComputeGEV(const std::vector<Eigen::Matrix2cd>& Cts, int t0)
{
    std::vector<GEVPPlateau> result;

    // for(int t = 0; t < t0; ++t) {
    // 	Eigen::ComplexEigenSolver<Eigen::Matrix2cd> eigsolv(Cts[t0].inverse() * Cts[t], false);
    // 	if (eigsolv.info() != Eigen::Success) abort();
    // 	result.push_back(GEVPPlateau(real(eigsolv.eigenvalues()(0)), real(eigsolv.eigenvalues()(1))));
    // 	//std::cout << eigsolv.eigenvalues() << std::endl;
    // }
    // // t0 slot = (0, 0)
    // result.push_back(GEVPPlateau(0., 0.));
    for(int t = 0; t < t0; ++t) {
	Eigen::ComplexEigenSolver<Eigen::Matrix2cd> eigsolv(Cts[t0].inverse() * Cts[t], false);
	if (eigsolv.info() != Eigen::Success) abort();
	result.push_back(GEVPPlateau(real(eigsolv.eigenvalues()(0)), real(eigsolv.eigenvalues()(1))));
	//std::cout << eigsolv.eigenvalues() << std::endl;
    }
    for(int t = t0; t < Cts.size()-t0; ++t) {
	Eigen::ComplexEigenSolver<Eigen::Matrix2cd> eigsolv(Cts[t0].inverse() * Cts[t], false);
	if (eigsolv.info() != Eigen::Success) abort();
	result.push_back(GEVPPlateau(real(eigsolv.eigenvalues()(1)), real(eigsolv.eigenvalues()(0))));
	//std::cout << eigsolv.eigenvalues() << std::endl;
    }
    for(int t = Cts.size()-t0; t < Cts.size(); ++t) {
	Eigen::ComplexEigenSolver<Eigen::Matrix2cd> eigsolv(Cts[t0].inverse() * Cts[t], false);
	if (eigsolv.info() != Eigen::Success) abort();
	result.push_back(GEVPPlateau(real(eigsolv.eigenvalues()(1)), real(eigsolv.eigenvalues()(0))));
	//std::cout << eigsolv.eigenvalues() << std::endl;
    }
    
    return result;
}

std::vector<GEVPPlateau> ComputeGEVPPlateau(const std::vector<std::vector<Eigen::Matrix2cd> >& rscorr, int t0, int tf)
{
    assert(t0 < rscorr[0].size());

    std::vector<GEVPPlateau> result;
    result.reserve(rscorr.size());
    
    typedef ResampledFitData<BootstrapResampler, double, double, Identity, Zero, Zero> RSFitData;
    
    PlateauFitModel* model = new PlateauFitModel();
    //ConstantPlateauFitModel* model = new ConstantPlateauFitModel();

    //int tsize = (tf > t0)? tf-1 : tf; 
    int tsize = 5;
    int tmin = 0;
    std::vector<std::vector<double>> e1(tsize-tmin, std::vector<double>(rscorr.size())), e2(tsize-tmin, std::vector<double>(rscorr.size()));
    std::vector<std::vector<double>> t(tsize-tmin, std::vector<double>(rscorr.size()));

    //std::vector<double> tmp1, tmp2, tmp3, tmp4;
    for(int n = 0; n < rscorr.size(); ++n) {
	std::vector<GEVPPlateau> gevplat = ComputeGEV(rscorr[n], t0);

	if(tf > t0 && t0 > tmin) {
	    for(int i = tmin ; i < t0 ; ++i) {
        int ii = i-tmin;
		e1[ii][n] = std::log(gevplat[i].E1) / ((double)(t0 - i));
		e2[ii][n] = std::log(gevplat[i].E2) / ((double)(t0 - i));
		t[ii][n] = (double)i;
	    }
	    for(int i = t0+1 ; i < tf ; ++i) {
        int ii = i-tmin;
		e1[ii-1][n] = std::log(gevplat[i].E1) / ((double)(t0 - i));
		e2[ii-1][n] = std::log(gevplat[i].E2) / ((double)(t0 - i));
		t[ii-1][n] = (double)i;
	    }
	}
	else {
	    for(int i = tmin ; i < tf ; ++i) {
        int ii = i-tmin;
		e1[ii][n] = std::log(gevplat[i].E1) / ((double)(t0 - i));
		e2[ii][n] = std::log(gevplat[i].E2) / ((double)(t0 - i));
		t[ii][n] = (double)i;
	    }
	}
    }

    RSFitData* fdata1 = new RSFitData(e1, t);
    RSFitData* fdata2 = new RSFitData(e2, t);
    
    for(unsigned int s = 0; s < fdata1->NSamples(); ++s) {
	//std::cout << "Sample " << s << " : " << std::endl;
	// Fit E1
	//fdata1->Value(s).DisablePoint(t0);
	FitResult<double, double> fit1 = Fitter<Chi2Base, Mn2MigradMinimizer>::Fit(&fdata1->Value(s), model);
	//std::cout << "Fit E1 : " << std::endl << fit1 << std::endl;
	// Fit E2
	// fdata2->Value(s).DisablePoint(t0);
	FitResult<double, double> fit2 = Fitter<Chi2Base, Mn2MigradMinimizer>::Fit(&fdata2->Value(s), model);
	//std::cout << "Fit E2 : " << std::endl << fit2 << std::endl;

	// Store results
	double E1 = fabs(fit1.FittedParameters().Value(0));
	double E2 = fabs(fit2.FittedParameters().Value(0));
	result.push_back(GEVPPlateau(E1, E2));
	// tmp1.push_back(fit1.FittedParameters().Value(0));
	// tmp2.push_back(fit1.FittedParameters().Value(1));
	// tmp3.push_back(fit2.FittedParameters().Value(0));
	// tmp4.push_back(fit2.FittedParameters().Value(1));
    }
    // std::cout << std::endl << "Fit 1 : " << BootstrapResampler<double>::Mean(tmp1) << std::endl << BootstrapResampler<double>::Mean(tmp2) << std::endl;
    // std::cout << std::endl << "Fit 2 : " << BootstrapResampler<double>::Mean(tmp3) << std::endl << BootstrapResampler<double>::Mean(tmp4) << std::endl;

    delete fdata1;
    delete fdata2;
    delete model;

    return result;
}

std::vector<Eigen::Matrix2cd> ReadGEVPCorrelator(const GEVPCorrelatorHeader& h, const std::string& file)
{
    return LQCDA::ParseFile<std::vector<Eigen::Matrix2cd> >(file, h);
}

void ReadGEVPSampleHelper(std::vector<std::vector<Eigen::Matrix2cd> >* sample, const GEVPCorrelatorHeader& h, const std::string& file)
{
    sample->push_back(LQCDA::ParseFile<std::vector<Eigen::Matrix2cd> >(file, h));
}

std::vector<std::vector<Eigen::Matrix2cd> > ReadGEVPSample(const GEVPCorrelatorHeader& h, const std::string& flist)
{
    using namespace std::placeholders;

    std::vector<std::vector<Eigen::Matrix2cd> >* buf = new std::vector<std::vector<Eigen::Matrix2cd> >();
    
    ManifestFileReader list(flist.c_str());
    list.ForEach(std::bind(ReadGEVPSampleHelper, buf, h, _1));

    std::vector<std::vector<Eigen::Matrix2cd> > result(*buf);
    delete buf;
    
    return result;
}

std::vector<std::vector<double> > ResampleLocalCorrelators(const std::vector<std::vector<double> >& sample, unsigned int Nboot, LQCDA::RandGen::rg_state state)
{
    RandGen rng;
    rng.setState(state);
    
    std::vector<std::vector<double> > result(Nboot);

    unsigned int NSample = sample.size();
    double NSample_inv = 1. / (double)NSample;
    unsigned int index;
    
    for(int n = 0; n < Nboot; ++n) {
	index = rng.getUniformInt(NSample);
	std::vector<double> buf(sample[index]);
	for(int i = 1; i < NSample; ++i) {
	    index = rng.getUniformInt(NSample);
	    buf = buf + sample[index];
	}
	buf = buf * NSample_inv;

	result[n] = buf;
    }

    return result;
}

std::vector<double> ReadLocalCorrelator(const GEVPCorrelatorHeader& h, const std::string& file)
{
    return LQCDA::ParseFile<std::vector<double> >(file, h);
}

void ReadLocalSampleHelper(std::vector<std::vector<double>>* sample, const GEVPCorrelatorHeader& h, const std::string& file)
{
    sample->push_back(LQCDA::ParseFile<std::vector<double>>(file, h));
}

std::vector<std::vector<double> > ReadLocalSample(const GEVPCorrelatorHeader& h, const std::string& flist)
{
    using namespace std::placeholders;

    std::vector<std::vector<double> >* buf = new std::vector<std::vector<double> >();
    
    ManifestFileReader list(flist.c_str());
    list.ForEach(std::bind(ReadLocalSampleHelper, buf, h, _1));

    std::vector<std::vector<double> > result(*buf);
    delete buf;

    return result;
}

std::vector<double> ComputeLocalPlateau(const std::vector<std::vector<double> >& rscorr, int tsize)
{
    int nboot = rscorr.size();

    std::vector<double> result;
    result.reserve(nboot);
    
    typedef ResampledFitData<BootstrapResampler, double, double, Identity, Zero, Zero> RSFitData;
    
    PlateauFitModel* model = new PlateauFitModel();
    //ConstantPlateauFitModel* model = new ConstantPlateauFitModel();

    int tmin = 1;tsize=15;

    std::vector<std::vector<double>> e(tsize-tmin, std::vector<double>(nboot));
    std::vector<std::vector<double>> t(tsize-tmin, std::vector<double>(nboot));

    for(int n = 0; n < nboot; ++n) {
	std::vector<double> corr = rscorr[n];
	for(int i = tmin ; i < tsize ; ++i) {
        int ii = i-tmin;
	    e[ii][n] = std::log(corr[i] / corr[i+1]);
	    t[ii][n] = (double)i;
	}
    }

    RSFitData* fdata = new RSFitData(e, t);

    for(unsigned int s = 0; s < fdata->NSamples(); ++s) {
	//std::cout << "Sample " << s << " : " << std::endl;
	// Fit E
	FitResult<double, double> fit = Fitter<Chi2Base, Mn2MigradMinimizer>::Fit(&fdata->Value(s), model);
	//std::cout << "Fit E1 : " << std::endl << fit << std::endl;

	// Store results
	double E = fabs(fit.FittedParameters().Value(0));
	result.push_back(E);
    }

    delete fdata;
    delete model;

    return result;
}

double DeltaM2_COM(double W, double MPi, double L)
{
    double E_COM = W;
    double q2 = (E_COM*E_COM/4. - MPi*MPi)*(L*L)/(4.*M_PI*M_PI);
    double delta;
    if(q2 < 0.) {
	delta = 0.;
	//std::cout << "q2 < 0 !" << std::endl;
    }
    else
	delta = -8.*M_PI*M_SQRTPI/(L*L*L)*q2/E_COM*Z001::z001q2(q2);
    // std::cout << "q2 = " << q2 << std::endl;
    return delta;
}

double DeltaM2_MV(double W, double MPi, double L, int d [3])
{
    double gamma = W / sqrt(fabs(W*W - 4.*M_PI*M_PI / (L*L) * (d[0]*d[0]+d[1]*d[1]+d[2]*d[2])));
    double E_COM = W / gamma;
    double q2 = (E_COM*E_COM/4. - MPi*MPi)*(L*L)/(4.*M_PI*M_PI);
    double delta;
    if(q2 < 0.)
	delta = 0.;
    else
	delta = -8.*M_PI*M_SQRTPI/(L*L*L)*q2/E_COM*Z001::z001q2(q2, gamma, d)/gamma;

    return delta;
}
