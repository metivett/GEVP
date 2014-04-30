/*
 * mass.cpp
 *
 *  Created on: Feb 18, 2013
 *      Author: Thibaut Metivet
 */

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include "LocalLuscherFit.hpp"
#include "GEVP.hpp"
#include "GEVPFitModels.hpp"
#include "Z001.hpp"
#include "LQCDA.hpp"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

 using namespace std;
 using namespace LQCDA;
 using namespace boost::filesystem;
 namespace po = boost::program_options;

 double GeV(double x, double beta)
 {
    // Map to convert lattice units to physical units
    map<double, double> ainv{{3.31, 1.697}, {3.5, 2.131}, {3.61, 2.561}, {3.7, 3.026}, {3.8, 3.662}};

    return x*ainv[beta];
}

void analyze(const char* fname, int t0, int nBoot, double L, int Nt, const char* tstr, double t1 = 0., double beta = 0.)
{
    // List-file for analyzed set
    std::string analyzed_file(fname);

    // Parameters
    enum Type { COM, MV };
    Type type;
    if(strcmp(tstr, "COM") == 0) type = Type::COM;
    else if(strcmp(tstr, "MV") == 0) type = Type::MV;
    else {
       cerr << "Type should be COM or MV !" << endl;
       return;
   }
   int d [3] = {1, 0, 0};
   int Tsize = Nt / 2;

    // Print info
   cout << endl << endl
   << "Performing GEVP analysis with :" << endl
   << "File " << analyzed_file << endl
   << "t0 = " << t0 << endl
   << "nBoot = " << nBoot << endl
   << "L = " << L << endl
   << "type " << tstr << endl
   << "t1 = " << t1 << endl << endl;
   
    // Generate random number generator state
   RandGen rng;
   RandGen::rg_state state;
   rng.getState(state);

   
   cout << "Results :" << endl << endl;
   
    /*** Local pion operator analysis (MV and COM) ***/

    // Time range
   int LOCPiTrg = (t1 > 0.)? t1 / 2 : Nt / 2;
   
    // Read PionPion correlators for all configurations
   GEVPCorrelatorHeader pi_h_COM("S", "S", "P", "P", 0, -1);
    //GEVPCorrelatorHeader pi_h_MV("S", "S", "P", "P", 1, -1);
   GEVPCorrelatorHeader pi_h_MV("S", "S", "P", "P", 0, -1);

   std::vector<std::vector<double> > LOCPionCorrAll_COM = ReadLocalSample(pi_h_COM, analyzed_file);
   std::vector<std::vector<double> > LOCPionCorrAll_MV = ReadLocalSample(pi_h_MV, analyzed_file);
   
    // Resample PP correlators (2000 bootstraps)
   std::vector<std::vector<double> > rsLOCPionCorr_COM = ResampleLocalCorrelators(LOCPionCorrAll_COM, nBoot, state);
   std::vector<std::vector<double> > rsLOCPionCorr_MV = ResampleLocalCorrelators(LOCPionCorrAll_MV, nBoot, state);

    // Effective mass
   std::vector<std::vector<double>> rsLOCPionMeff_COM(nBoot);
   for(int s = 0; s < nBoot; ++s) {
       rsLOCPionMeff_COM[s].reserve(Tsize);
       for(int i = 0; i < Nt; ++i) {
           rsLOCPionMeff_COM[s].push_back(std::log(rsLOCPionCorr_COM[s][i] / rsLOCPionCorr_COM[s][i+1]));
           if(!std::isfinite(rsLOCPionMeff_COM[s][i]))
              LOCPiTrg = std::min(i - 1, LOCPiTrg);
      };
  }
  std::vector<std::vector<double>> rsLOCPionMeff_MV(nBoot);
  for(int s = 0; s < nBoot; ++s) {
   rsLOCPionMeff_MV[s].reserve(Tsize);
   for(int i = 0; i < Nt; ++i) {
       rsLOCPionMeff_MV[s].push_back(std::log(rsLOCPionCorr_MV[s][i] / rsLOCPionCorr_MV[s][i+1]));
       if(!std::isfinite(rsLOCPionMeff_MV[s][i]))
          LOCPiTrg = std::min(i - 1, LOCPiTrg);
  };
}

std::vector<double> LOCPionMeff_COM = BootstrapResampler::Mean(rsLOCPionMeff_COM);
std::vector<double> LOCPionMefferr_COM = BootstrapResampler::Variance(rsLOCPionMeff_COM);
std::vector<double> LOCPionMeff_MV = BootstrapResampler::Mean(rsLOCPionMeff_MV);
std::vector<double> LOCPionMefferr_MV = BootstrapResampler::Variance(rsLOCPionMeff_MV);

    // Print effective mass
ofstream LOCPionMeffof_COM("LOCPionMeff_COM.dat");
for(int i = 0; i < LOCPionMeff_COM.size(); ++i) {
	LOCPionMeffof_COM << i << ' ' << LOCPionMeff_COM[i] << ' ' << sqrt(LOCPionMefferr_COM[i]) << endl;
}
ofstream LOCPionMeffof_MV("LOCPionMeff_MV.dat");
for(int i = 0; i < LOCPionMeff_MV.size(); ++i) {
	LOCPionMeffof_MV << i << ' ' << LOCPionMeff_MV[i] << ' ' << sqrt(LOCPionMefferr_MV[i]) << endl;
}

    // Compute resampled plateaus
std::vector<double> rsLOCPionPlat_COM = ComputeLocalPlateau(rsLOCPionCorr_COM, LOCPiTrg);
std::vector<double> rsLOCPionPlat_MV = ComputeLocalPlateau(rsLOCPionCorr_MV, LOCPiTrg);

    // std::vector<double> rsLOCPionPlat_COM(nBoot);
    // for(int s = 0; s < nBoot; ++s) {
    // 	rsLOCPionPlat_COM[s] = sqrt(rsLOCPionPlat[s]*rsLOCPionPlat[s] - 4.*M_PI*M_PI/(L*L)*(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]));
    // }

double LOCPionPlat_COM = BootstrapResampler::Mean(rsLOCPionPlat_COM);
double LOCPionPlaterr_COM = sqrt(BootstrapResampler::Variance(rsLOCPionPlat_COM));
double LOCPionPlat_MV = BootstrapResampler::Mean(rsLOCPionPlat_MV);
double LOCPionPlaterr_MV = sqrt(BootstrapResampler::Variance(rsLOCPionPlat_MV));

    // Print results
cout << "Pion analysis : " << endl
<< "tf = " << LOCPiTrg << endl
<< "E_MV = " << LOCPionPlat_MV << " (" << GeV(LOCPionPlat_MV, beta) << " GeV)"
<< " +- " << LOCPionPlaterr_MV << " (" << GeV(LOCPionPlaterr_MV, beta) << " GeV)"
<< endl
<< "E_COM = " << LOCPionPlat_COM << " (" << GeV(LOCPionPlat_COM, beta) << " GeV)"
<< " +- " << LOCPionPlaterr_COM << " (" << GeV(LOCPionPlaterr_COM, beta) << " GeV)"
<< endl << endl;

cout << endl;



    /*** GEVP analysis ***/

    // Time range
int GEVPTrg = 2*t0;

    // Read GEVP correlator matrices for all configurations
GEVPCorrelatorHeader GEVP_h("S", "S", "Vi", "Vi", 1, -1);
std::vector<std::vector<Eigen::Matrix2cd> > GEVPcorrAll = ReadGEVPSample(GEVP_h, analyzed_file);

    // Resample GEVP correlator matrices (2000 boostraps)
std::vector<std::vector<Eigen::Matrix2cd> > rsGEVPcorr = ResampleGEVPCorrelators(GEVPcorrAll, nBoot, state);

    // //Treat t1 trick
if(t1 > 0.) {
   for(int s = 0; s < rsGEVPcorr.size(); ++s) {
       for(int t = 0; t < rsGEVPcorr[s].size(); ++t) {
          double tmp = exp(rsLOCPionPlat_MV[s] * (t1-(double)t));
          rsGEVPcorr[s][t](0,1) *= tmp;
          rsGEVPcorr[s][t](1,1) *= tmp;
      }
  }
}

    // Print mean(rsGEVPcorr)
std::vector<Eigen::Matrix2cd> GEVPCorr = BootstrapResampler::Mean(rsGEVPcorr);
ofstream GEVPCorrof("GEVPCorr.dat");
for(int i = 0; i < GEVPCorr.size(); ++i) {
	GEVPCorrof << i << ' ' << GEVPCorr[i] << endl;
}

    // Generalized eigenvalues
std::vector<std::vector<GEVPPlateau>> rsGEV(rsGEVPcorr.size());
for(int s = 0; s < rsGEV.size(); ++s) {
	rsGEV[s] = ComputeGEV(rsGEVPcorr[s], t0);
}
    // Mean and variance
std::vector<GEVPPlateau> GEV = BootstrapResampler::Mean(rsGEV);
std::vector<GEVPPlateau> GEVerr = BootstrapResampler::Variance(rsGEV);
    // Print central value and error
ofstream gevof1("GEV1.dat");
ofstream gevof2("GEV2.dat");
for(int i = 0; i < GEV.size(); ++i) {
	gevof1 << i << ' ' << GEV[i].E1 << ' ' << sqrt(GEVerr[i].E1) << endl;
	gevof2 << i << ' ' << GEV[i].E2 << ' ' << sqrt(GEVerr[i].E2) << endl;
}

    // Effective masses
std::vector<std::vector<GEVPPlateau>> rsGEVPMeff(rsGEVPcorr.size());
for(int s = 0; s < nBoot; ++s) {
	for(int i = 0; i < t0; ++i) {
       rsGEVPMeff[s].push_back(GEVPPlateau(std::log(rsGEV[s][i].E1)/(double)(t0-i), std::log(rsGEV[s][i].E2)/(double)(t0-i)));
       if(!std::isfinite(rsGEVPMeff[s][i].E1) || !std::isfinite(rsGEVPMeff[s][i].E2))
          GEVPTrg = std::min(i - 1, GEVPTrg);
  }
  rsGEVPMeff[s].push_back(GEVPPlateau(0., 0.));
  for(int j = t0+1; j < 2*t0; ++j) {
   rsGEVPMeff[s].push_back(GEVPPlateau(std::log(rsGEV[s][j].E1)/(double)(t0-j), std::log(rsGEV[s][j].E2)/(double)(t0-j)));
   if(!std::isfinite(rsGEVPMeff[s][j].E1) || !std::isfinite(rsGEVPMeff[s][j].E2))
      GEVPTrg = std::min(j - 1, GEVPTrg);
}
}
    // Mean and variance
std::vector<GEVPPlateau> GEVPMeff = BootstrapResampler::Mean(rsGEVPMeff);
std::vector<GEVPPlateau> GEVPMefferr = BootstrapResampler::Variance(rsGEVPMeff);
    // Print central value and error
ofstream GEVPMeffof1("GEVPMeff1.dat");
ofstream GEVPMeffof2("GEVPMeff2.dat");
for(int i = 0; i < GEVPMeff.size(); ++i) {
	GEVPMeffof1 << i << ' ' << GEVPMeff[i].E1 << ' '
  << sqrt(GEVPMefferr[i].E1) << endl;
  GEVPMeffof2 << i << ' ' << GEVPMeff[i].E2 << ' '
  << sqrt(GEVPMefferr[i].E2) << endl;
}

    // ofstream matof("mat2.dat");
    // for(int i = 0; i < rsGEVPcorr[0].size(); ++i) {
    // 	matof << i << ' '
    // 	      << real(rsGEVPcorr[0][i](0,0)) << ' '
    // 	      << imag(rsGEVPcorr[0][i](0,1)) << ' '
    // 	      << imag(rsGEVPcorr[0][i](1,0)) << ' '
    // 	      << real(rsGEVPcorr[0][i](1,1)) << endl;
    // } 

    // Resampled Plateaus
if(GEVPTrg > 2*t0) GEVPTrg = 2*t0;
std::vector<GEVPPlateau> rsGEVPPlat = ComputeGEVPPlateau(rsGEVPcorr, t0, GEVPTrg);

GEVPPlateau GEVPPlat = BootstrapResampler::Mean(rsGEVPPlat);
GEVPPlateau GEVPPlaterr = sqrt(BootstrapResampler::Variance(rsGEVPPlat));

cout << "GEVP analysis : " << endl
<< "tf = " << GEVPTrg << endl
<< "Rho1 #0 : " << GEVPPlat.E1 << " (" << GeV(GEVPPlat.E1, beta) << " GeV)"
<< " +- " << GEVPPlaterr.E1 << " (" << GeV(GEVPPlaterr.E1, beta) << " GeV)"
<< endl
<< "Rho1 #1 : " << GEVPPlat.E2 << " (" << GeV(GEVPPlat.E2, beta) << " GeV)"
<< " +- " << GEVPPlaterr.E2 << " (" << GeV(GEVPPlaterr.E2, beta) << " GeV)"
<< endl << endl;



    /*** Local rho operator analysis ***/

    // Time range
int LOCRhoTrg = (t1 > 0.)? t1 / 2 : Nt / 2;

    // Get resampled rho-rho correlators (2000 bootstraps)
std::vector<std::vector<double> > rsLOCRhoCorr_MV(nBoot);
for(int n = 0; n < nBoot; ++n) {
	rsLOCRhoCorr_MV.reserve(Nt);
	for(int i = 0; i < Nt; ++i) {
       rsLOCRhoCorr_MV[n].push_back(std::real(rsGEVPcorr[n][i](0,0)));
   };
}
GEVPCorrelatorHeader rho_h_COM("S", "S", "Vi", "Vi", 0, -1);
std::vector<std::vector<double> > LOCRhoCorrAll_COM = ReadLocalSample(rho_h_COM, analyzed_file);
std::vector<std::vector<double> > rsLOCRhoCorr_COM = ResampleLocalCorrelators(LOCRhoCorrAll_COM, nBoot, state);

    // Effective mass
std::vector<std::vector<double>> rsLOCRhoMeff_COM(nBoot);
std::vector<std::vector<double>> rsLOCRhoMeff_MV(nBoot);
for(int s = 0; s < nBoot; ++s) {
	rsLOCRhoMeff_COM[s].reserve(Tsize);
	for(int i = 0; i < Nt; ++i) {
       rsLOCRhoMeff_COM[s].push_back(std::log(rsLOCRhoCorr_COM[s][i] / rsLOCRhoCorr_COM[s][i+1]));
       if(!std::isfinite(rsLOCRhoMeff_COM[s][i]))
          LOCRhoTrg = std::min(i - 1, LOCRhoTrg);
  };
}
for(int s = 0; s < nBoot; ++s) {
	rsLOCRhoMeff_MV[s].reserve(Tsize);
	for(int i = 0; i < Nt; ++i) {
       rsLOCRhoMeff_MV[s].push_back(std::log(rsLOCRhoCorr_MV[s][i] / rsLOCRhoCorr_MV[s][i+1]));
       if(std::isnan(rsLOCRhoMeff_MV[s][i]))
          LOCRhoTrg = std::min(i - 1, LOCRhoTrg);
  };
}

std::vector<double> LOCRhoMeff_COM = BootstrapResampler::Mean(rsLOCRhoMeff_COM);
std::vector<double> LOCRhoMefferr_COM = BootstrapResampler::Variance(rsLOCRhoMeff_COM);

std::vector<double> LOCRhoMeff_MV = BootstrapResampler::Mean(rsLOCRhoMeff_MV);
std::vector<double> LOCRhoMefferr_MV = BootstrapResampler::Variance(rsLOCRhoMeff_MV);

    // Print effective mass
ofstream LOCRhoMeffof_COM("LOCRhoMeff_COM.dat");
for(int i = 0; i < LOCRhoMeff_COM.size(); ++i) {
	LOCRhoMeffof_COM << i << ' ' << LOCRhoMeff_COM[i] << ' ' << sqrt(LOCRhoMefferr_COM[i]) << endl;
}
ofstream LOCRhoMeffof_MV("LOCRhoMeff_MV.dat");
for(int i = 0; i < LOCRhoMeff_MV.size(); ++i) {
	LOCRhoMeffof_MV << i << ' ' << LOCRhoMeff_MV[i] << ' ' << sqrt(LOCRhoMefferr_MV[i]) << endl;
}

    // Compute resampled plateaus
std::vector<double> rsLOCRhoPlat_COM = ComputeLocalPlateau(rsLOCRhoCorr_COM, LOCRhoTrg);
std::vector<double> rsLOCRhoPlat_MV = ComputeLocalPlateau(rsLOCRhoCorr_MV, LOCRhoTrg);

double LOCRhoPlat_COM = BootstrapResampler::Mean(rsLOCRhoPlat_COM);
double LOCRhoPlaterr_COM = sqrt(BootstrapResampler::Variance(rsLOCRhoPlat_COM));
double LOCRhoPlat_MV = BootstrapResampler::Mean(rsLOCRhoPlat_MV);
double LOCRhoPlaterr_MV = sqrt(BootstrapResampler::Variance(rsLOCRhoPlat_MV));

    // Print results
cout << "Rho analysis : " << endl
<< "tf = " << LOCRhoTrg << endl
<< "E_MV = " << LOCRhoPlat_MV << " (" << GeV(LOCRhoPlat_MV, beta) << " GeV)"
<< " +- " << LOCRhoPlaterr_MV << " (" << GeV(LOCRhoPlaterr_MV, beta) << " GeV)"
<< endl
<< "E_COM = " << LOCRhoPlat_COM << " (" << GeV(LOCRhoPlat_COM, beta) << " GeV)"
<< " +- " << LOCRhoPlaterr_COM << " (" << GeV(LOCRhoPlaterr_COM, beta) << " GeV)"
<< endl << endl;



    // Luscher analysis
std::vector<double> rsMRho2(nBoot);
std::vector<double> rsg2(nBoot);
std::vector<double> rsMRho(nBoot);
std::vector<double> rsg(nBoot);

switch(type) {
    case Type::MV:
    for(int s = 0; s < nBoot; ++s) {
       double E2_1_com = rsGEVPPlat[s].E1*rsGEVPPlat[s].E1 - 4.*M_PI*M_PI / (L*L) * (d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
       double E2_2_com = rsGEVPPlat[s].E2*rsGEVPPlat[s].E2 - 4.*M_PI*M_PI / (L*L) * (d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
       double MPicom = rsLOCPionPlat_COM[s];
       double DeltaM2_1 = DeltaM2_MV(rsGEVPPlat[s].E1, MPicom, L, d);
       double DeltaM2_2 = DeltaM2_MV(rsGEVPPlat[s].E2, MPicom, L, d);

       double MRho2 = fabs((E2_1_com*DeltaM2_2 - E2_2_com*DeltaM2_1)/(DeltaM2_2 - DeltaM2_1));
       double g2 = 6.*M_PI*fabs((E2_1_com-E2_2_com)/(DeltaM2_1-DeltaM2_2));

       rsMRho2[s] = MRho2;
       rsg2[s] = g2;

       rsMRho[s] = sqrt(MRho2);
       rsg[s] = sqrt(g2);
	    // cout << "gamma = " << gamma << endl;
	    // cout << "MPicom = " << MPicom << endl;
	    // cout << "g2 = " << g2 << endl;
	    // cout << "Mrho2 = " << MRho2 << endl;
   }
   break;
   
   case Type::COM:
   for(int s = 0; s < nBoot; ++s) {
       double E2_1_com = rsGEVPPlat[s].E1*rsGEVPPlat[s].E1;
       double E2_2_com = rsGEVPPlat[s].E2*rsGEVPPlat[s].E2;
       double MPicom = rsLOCPionPlat_COM[s];
       double DeltaM2_1 = DeltaM2_COM(rsGEVPPlat[s].E1, MPicom, L);
       double DeltaM2_2 = DeltaM2_COM(rsGEVPPlat[s].E2, MPicom, L);

       double MRho2 = fabs((E2_1_com*DeltaM2_2 - E2_2_com*DeltaM2_1)/(DeltaM2_2 - DeltaM2_1));
       double g2 = 6.*M_PI*fabs((E2_1_com-E2_2_com)/(DeltaM2_1-DeltaM2_2));

       rsMRho2[s] = MRho2;
       rsg2[s] = g2;

       rsMRho[s] = sqrt(MRho2);
       rsg[s] = sqrt(g2);

	    // cout << "g2 = " << g2 << endl;
	    // cout << "Mrho2 = " << MRho2 << endl;
   }
   break;
}

double g = BootstrapResampler::Mean(rsg);
double gerr = sqrt(BootstrapResampler::Variance(rsg));
double MRho = BootstrapResampler::Mean(rsMRho);
double MRhoerr = sqrt(BootstrapResampler::Variance(rsMRho));

cout << "Coupling parameters : " << endl
<< "g2 = " << BootstrapResampler::Mean(rsg2)
<< "  +- " <<  sqrt(BootstrapResampler::Variance(rsg2))
<< endl
<< "Mrho2 = " << BootstrapResampler::Mean(rsMRho2)
<< "  +- " <<  sqrt(BootstrapResampler::Variance(rsMRho2))
<< endl << endl
<< "g = " << g
<< " +- " <<  gerr
<< endl
<< "Mrho = " << MRho << " (" << GeV(MRho, beta) << " GeV)"
<< " +- " << MRhoerr << " (" << GeV(MRhoerr, beta) << " GeV)"
<< endl << endl;

    // Print summary into file and resampled data in other file
string sumfstr;
string rsdatfstr;

path p(analyzed_file);
if(p.has_parent_path()) {
	path parent = p.parent_path();
	if(parent.has_parent_path()) {
       sumfstr = parent.parent_path().filename().native() + "-" + parent.filename().native();
       rsdatfstr = sumfstr + ".rsdat";
   }
   else {
       sumfstr = string("summary.dat");
       rsdatfstr = string("summary.rsdat");
   }
}
else {
	sumfstr = string("summary.dat");
	rsdatfstr = string("summary.rsdat");
}

ofstream sumof(sumfstr);
sumof << "Performing GEVP analysis with :" << endl
<< "File " << analyzed_file << endl
<< "t0 = " << t0 << endl
<< "nBoot = " << nBoot << endl
<< "L = " << L << endl
<< "type " << tstr << endl << endl << endl;

sumof << "Results :" << endl << endl;

sumof << "GEVP analysis : " << endl
<< "tf = " << GEVPTrg << endl
<< "Rho1 #0 : " << GEVPPlat.E1 << " (" << GeV(GEVPPlat.E1, beta) << " GeV)"
<< " +- " << GEVPPlaterr.E1 << " (" << GeV(GEVPPlaterr.E1, beta) << " GeV)"
<< endl
<< "Rho1 #1 : " << GEVPPlat.E2 << " (" << GeV(GEVPPlat.E2, beta) << " GeV)"
<< " +- " << GEVPPlaterr.E2 << " (" << GeV(GEVPPlaterr.E2, beta) << " GeV)"
<< endl << endl;
sumof << "Rho analysis : " << endl
<< "tf = " << LOCRhoTrg << endl
<< "E_MV = " << LOCRhoPlat_MV << " (" << GeV(LOCRhoPlat_MV, beta) << " GeV)"
<< " +- " << LOCRhoPlaterr_MV << " (" << GeV(LOCRhoPlaterr_MV, beta) << " GeV)"
<< endl
<< "E_COM = " << LOCRhoPlat_COM << " (" << GeV(LOCRhoPlat_COM, beta) << " GeV)"
<< " +- " << LOCRhoPlaterr_COM << " (" << GeV(LOCRhoPlaterr_COM, beta) << " GeV)"
<< endl << endl;
sumof << "Pion analysis : " << endl
<< "tf = " << LOCRhoTrg << endl
<< "E_MV = " << LOCPionPlat_MV << " (" << GeV(LOCPionPlat_MV, beta) << " GeV)"
<< " +- " << LOCPionPlaterr_MV << " (" << GeV(LOCPionPlaterr_MV, beta) << " GeV)"
<< endl
<< "E_COM = " << LOCPionPlat_COM << " (" << GeV(LOCPionPlat_COM, beta) << " GeV)"
<< " +- " << LOCPionPlaterr_COM << " (" << GeV(LOCPionPlaterr_COM, beta) << " GeV)"
<< endl << endl;
sumof << "Coupling parameters : " << endl
<< "g = " << g
<< " +- " <<  gerr
<< endl
<< "Mrho = " << MRho << " (" << GeV(MRho, beta) << " GeV)"
<< " +- " << MRhoerr << " (" << GeV(MRhoerr, beta) << " GeV)"
<< endl << endl;

ofstream rssumof(rsdatfstr);
rssumof << "MpiCOM" << setw(15) << "MpiMV" << setw(15)
<< "RhoGEVP #0" << setw(15) << "RhoGEVP #1" << setw(15)
<< "RhoCOM" << setw(15) << "RhoMV" << setw(15)
<< "g2" << setw(15) << "g" << setw(15)
<< "Mrho2" << setw(15) << "Mrho"
<< endl;
for(int s = 0; s < nBoot; ++s) {
	rssumof << rsLOCPionPlat_COM[s] << setw(15) << rsLOCPionPlat_MV[s] << setw(15)
  << rsGEVPPlat[s].E1 << setw(15) << rsGEVPPlat[s].E2 << setw(15)
  << rsLOCRhoPlat_COM[s] << setw(15) << rsLOCRhoPlat_MV[s] << setw(15)
  << rsg2[s] << setw(15) << rsg[s] << setw(15)
  << rsMRho2[s] << setw(15) << rsMRho[s]
  << endl;
}
}

int main(int argc, char **argv) {
    po::positional_options_description p;
    p.add("input-file", -1);
    
    po::options_description opts("Usage");
    opts.add_options()
    ("help", "print help message")
    ("inputlist", po::value<std::string>(), "process list of input files")
    ("inputfile", po::value<std::string>(), "process a single file")
    ("t0", po::value<int>(), "set the t0 reference time")
    ("L", po::value<int>(), "space extent of the lattice")
    ("Nt", po::value<int>(), "time extent of the lattice")
    ("type", po::value<string>(), "set type of analysis (COM or MV)")
    ("t1", po::value<int>(), "if t1 > 0, apply t1-trick exponential correction")
    ("beta", po::value<double>(), "provide optional beta value to convert lattices units into physical units")
    ("nboot", po::value<int>(), "set number of bootstraps used for statistical resampling analysis")
    ;
    
    int nBoot = 2000;
    
    // Get arguments
    if(argc == 7) {
       analyze(argv[1], atoi(argv[2]), nBoot, atof(argv[3]), atoi(argv[4]), argv[5], atof(argv[6]));
       return 0;
   }
   else if(argc == 8) {
       analyze(argv[1], atoi(argv[2]), nBoot, atof(argv[3]), atoi(argv[4]), argv[5], atof(argv[6]), atof(argv[7]));
       return 0;
   }
   else if(argc == 2) {
       ifstream manfile(argv[1]);
       string buf;
       while(getline(manfile, buf)) {
           istringstream input(buf);
           string fname;
           int t0;
           double L;
           int Nt;
           string tstr;
           double t1;
           input >> fname;
           input >> t0;
           input >> L;
           input >> Nt;
           input >> tstr;
           input >> t1;
           if(input) {
              double beta;
              input >> beta;
              analyze(fname.c_str(), t0, nBoot, L, Nt, tstr.c_str(), t1, beta);
          }
          else {
              analyze(fname.c_str(), t0, nBoot, L, Nt, tstr.c_str(), t1);
          }
      }
      return 0;
  }
  else {
   cout << "Usage : LQCDUtils filename t0 L Nt type(=MV or COM) t1 [beta]" << endl
   << "     or LQCDUtils filename [beta]" << endl;
   return 1;
}

}
