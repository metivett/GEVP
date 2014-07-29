/*
 * Z001.cpp
 *
 *  Created on: Ju 05, 2013
 *      Author: Thibaut Metivet
 */

#include "Z001.hpp"
#include <gsl/gsl_integration.h>
#include <algorithm>

using namespace std;

namespace {

static const int N2MAX = 2501;	// was 10001
static const double EPSILON = 1.e-10;

class Array {
private:
    int * data;

    typedef int (*InitFcn)(unsigned int);
public:
    Array(unsigned int n, InitFcn f) :
	data(new int[n])
	{
	    for(unsigned int i = 0; i < n; ++i)
		data[i] = f(i);
	}

    int& operator[](unsigned int i) { return data[i]; }
};

bool is_integer(double x)
{
    if(fmod(x,floor(x))==0.)
	return true;
    else
	return false;
}

int init_n2degen(unsigned int nsq) {
    int n1,n2;
    int nsqmax,nsqmaxovrt2,sum;
    double rtnsq,rtdum;

    sum=0;
    rtnsq=sqrt(nsq);
    nsqmax=floor(rtnsq);
    nsqmaxovrt2=floor(sqrt(nsq/2.0));
    if (nsq==0) sum=1;
    else
    {
        if (is_integer(rtnsq)) sum=6;
        else sum=0;
        for(n1=1;n1<=nsqmax;n1++)
        {
            rtdum=sqrt(nsq-n1*n1);
            if(is_integer(rtdum) && rtdum !=0.) sum+=12;
        }
        for(n1=1;n1<=nsqmaxovrt2;n1++)
        {
            rtdum=sqrt(nsq-2.0*n1*n1);
            if(is_integer(rtdum) && rtdum !=0.) sum+=8;
        }
        for(n1=1;n1<=nsqmax;n1++)
        {
            for(n2=1;n2<=min(n1-1.0,floor(sqrt(nsq-n1*n1)));n2++)
            {
                rtdum=sqrt(nsq-n1*n1-n2*n2);
                if(is_integer(rtdum) && rtdum !=0.) sum+=16;
            }
        }
    }
    return sum;
}

double fnc1(double t, void * params)
{
    double q2 = *(double*) params;
    return (exp(q2*t)-1.0)/pow(t,1.5);
}

double fnc2(double t, void * params)
{
    double * par = (double*) params;
    double q2 = *par;
    double n2 = *(++par);
    return exp(q2*t-M_PI*M_PI*n2/t)/pow(t,1.5);	
}


double integral(const gsl_function* F, double a, double b)
{
    gsl_integration_workspace * w 
	= gsl_integration_workspace_alloc (1000);
       
       double result, error;

       double eps = 1.e-7;
     
       gsl_integration_qags (F, a, b, eps, eps, 1000,
                             w, &result, &error); 
     
       gsl_integration_workspace_free (w);
       
       return result;
}

double sum1(double qsq, int N2MAX, double EPS)
{
    double sum, dum;
    int nsq;

    sum = 0.;
    nsq = 0;
    do
    {
	dum = Z001::n2degen(nsq)*expf(qsq-nsq)/(nsq-qsq);
	sum += dum;
    } while (++nsq < N2MAX && (fabsf(dum) > fabsf(EPS*sum) || dum == 0.));
    if (nsq == N2MAX) printf("sum1 did not converge!\n");

    return sum;
}

void GetNParallel(double * npar, int n[3], int d[3])
{
    double d2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
    if(d2 == 0.) {
	npar[0] = 0.;
	npar[1] = 0.;
	npar[2] = 0.;
    }
    else {
	double ndotd = n[0]*d[0] + n[1]*d[1] + n[2]*d[2];
	double ndotd_d2 = ndotd / d2;
	npar[0] = ndotd_d2 * (double)d[0];
	npar[1] = ndotd_d2 * (double)d[1];
	npar[2] = ndotd_d2 * (double)d[2];
    }
}

double sum1(double qsq, double gamma, int d [3], int N2MAX, double EPS)
{
    double sum, dum, r2;

    int n [3];
    double n_parallel [3], n_ortho[3], buf[3];

    int nmax= sqrt(N2MAX);

    sum = 0.;
    for(int n1 = -nmax; n1 < nmax; ++n1)
	for(int n2 = -sqrt(N2MAX-n1*n1); n2 < sqrt(N2MAX-n1*n1); ++n2)
	    for(int n3 = -sqrt(N2MAX-n1*n1-n2*n2); n3 < sqrt(N2MAX-n1*n1-n2*n2); ++n3) {
		n[0] = n1; n[1] = n2; n[2] = n3;
		GetNParallel(n_parallel, n, d);
		n_ortho[0] = n[0] - n_parallel[0];
		n_ortho[1] = n[1] - n_parallel[1];
		n_ortho[2] = n[2] - n_parallel[2];
		double n_ortho2 = n_ortho[0] * n_ortho[0] + n_ortho[1] * n_ortho[1] + n_ortho[2] * n_ortho[2];
		buf[0] = (double)n_parallel[0] + (double)d[0] / 2.;
		buf[1] = (double)n_parallel[1] + (double)d[1] / 2.;
		buf[2] = (double)n_parallel[2] + (double)d[2] / 2.;
		double buf2 = buf[0] * buf[0] + buf[1] * buf[1] + buf[2] * buf[2];
		r2 = buf2 / (gamma * gamma) + n_ortho2;
		
		dum = exp(qsq-r2) / (r2 - qsq);
		sum += dum;
	    }
    if(fabsf(dum) < fabsf(EPS*sum))
        return sum;

    printf("sum1 did not converge!\n");

    return sum;
}

double sum2(double qsq, int N2MAX, double EPS)
{
  double sum,dum;

  gsl_function F2;
  F2.function = &fnc2;
  
  double par[2];
  par[0] = qsq;
  F2.params = &par;

  sum=0.;
  int nsq=1;
  do
  {
      par[1] = (double)nsq;
      dum = Z001::n2degen(nsq)*integral(&F2,0.0,1.0);
      sum += dum;
  } while (++nsq < N2MAX && (fabs(dum) > fabs(EPS*sum) || dum == 0.));
  if (nsq == N2MAX) printf("sum2 did not converge!\n");
  return sum;
}

double sum2(double qsq, double gamma, int d [3], int N2MAX, double EPS)
{
    gsl_function F2;
    F2.function = &fnc2;
    
    double par[2];
    par[0] = qsq;
    F2.params = &par;

    double sum, dum, gn2;

    int n [3];
    double n_parallel [3], n_ortho[3], buf[3];

    int nmax= sqrt(N2MAX);

    sum = 0.;
    for(int n1 = -nmax; n1 < nmax; ++n1)
	for(int n2 = -sqrt(N2MAX-n1*n1); n2 < sqrt(N2MAX-n1*n1); ++n2)
	    for(int n3 = -sqrt(N2MAX-n1*n1-n2*n2); n3 < sqrt(N2MAX-n1*n1-n2*n2); ++n3) {
		n[0] = n1; n[1] = n2; n[2] = n3;
		GetNParallel(n_parallel, n, d);
		double n_parallel2 = n_parallel[0] * n_parallel[0] + n_parallel[1] * n_parallel[1] + n_parallel[2] * n_parallel[2];
		n_ortho[0] = n[0] - n_parallel[0];
		n_ortho[1] = n[1] - n_parallel[1];
		n_ortho[2] = n[2] - n_parallel[2];
		double n_ortho2 = n_ortho[0] * n_ortho[0] + n_ortho[1] * n_ortho[1] + n_ortho[2] * n_ortho[2];
		gn2 = (gamma * gamma) * n_parallel2 + n_ortho2;
		if(gn2 == 0.) continue;
		par[1] = gn2;
		int ndotd = n[0]*d[0] + n[1]*d[1] + n[2]*d[2];

		if(ndotd % 2 == 0)
		    dum = integral(&F2,0.0,1.0);
		else
		    dum = - integral(&F2,0.0,1.0);
		
		sum += dum;	
	    }
    if(fabsf(dum) < fabsf(EPS*sum))
	return sum;
    
    printf("sum2 did not converge!\n");

    return sum;
}

}

int Z001::n2degen(unsigned int n2)
{
    static Array * _n2degen = new Array(N2MAX, &init_n2degen);

    if(n2 < N2MAX)
	return (*_n2degen)[n2];
    else
	return init_n2degen(n2);
}

double Z001::z001q2(double q2)
{    
    gsl_function F1;
    F1.function = &fnc1;
    F1.params = &q2;
    
    double res = - M_PI;
    res += 0.25 * M_2_SQRTPI * sum1(q2, N2MAX, EPSILON);
    res += M_PI_2 * integral(&F1,0.0,1.0);
    res += M_PI_2 * sum2(q2, N2MAX, EPSILON);

    return res;
}



double Z001::z001q2(double q2, double gamma, int direction [3])
{
    gsl_function F1;
    F1.function = &fnc1;
    F1.params = &q2;
    
    double res = - gamma * M_PI;
    res += 0.25 * M_2_SQRTPI * sum1(q2, gamma, direction, N2MAX, EPSILON);
    res += gamma * M_PI_2 * integral(&F1,0.0,1.0);
    res += gamma * M_PI_2 * sum2(q2, gamma, direction, N2MAX/10, EPSILON);

    return res;
}
