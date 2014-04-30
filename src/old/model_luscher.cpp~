/*
 * model_luscher.cpp
 *
 *  Created on: Apr 18, 2013
 *      Author: Thibaut Metivet
 */

#include "model_luscher.hpp"
#include <gsl/gsl_integration.h>

using namespace LocalLuscher;
using namespace LQCDA;
using namespace std;

const double LuscherModel::EPSILON = 1.e-10;
int LuscherModel::_n2degen[N2MAX];

bool is_integer(double x)
{
    if(fmod(x,floor(x))==0.)
	return true;
    else
	return false;
}

int LuscherModel::degen(size_t nsq)
{
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

void LuscherModel::build_n2degen()
{
    for(int nsq=0; nsq<N2MAX; nsq++)
    {
	LuscherModel::_n2degen[nsq] = degen(nsq);
    }
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

double sum1(double qsq, int nuofnsq[], int N2MAX, double EPS)
{
    double sum, dum;
    int nsq;

    sum = 0.;
    nsq = 0;
    do
    {
	dum = nuofnsq[nsq]*expf(qsq-nsq)/(nsq-qsq);
	sum += dum;
    } while (++nsq < N2MAX && (fabsf(dum) > fabsf(EPS*sum) || dum == 0.));
    if (nsq == N2MAX) printf("sum1 did not converge!\n");

    return sum;
}

double sum2(double qsq, int nuofnsq[], int N2MAX, double EPS)
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
      dum = nuofnsq[nsq]*integral(&F2,0.0,1.0);
      sum += dum;
  } while (++nsq < N2MAX && (fabs(dum) > fabs(EPS*sum) || dum == 0.));
  if (nsq == N2MAX) printf("sum2 did not converge!\n");
  return sum;
}



double LuscherModel::z001q2(double q2)
{
    gsl_function F1;
    F1.function = &fnc1;
    F1.params = &q2;
    
    double res = - M_PI;
    res += 0.25 * M_2_SQRTPI * sum1(q2, _n2degen, N2MAX, EPSILON);
    res += M_PI_2 * integral(&F1,0.0,1.0);
    res += M_PI_2 * sum2(q2, _n2degen, N2MAX, EPSILON);

    return res;
}
