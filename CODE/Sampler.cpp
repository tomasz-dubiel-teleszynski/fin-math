//Sampler.cpp
#include "Sampler.h"
#include<cstdlib>
#include<cmath>


//Uniform
UniformSampler::UniformSampler(double lower, double upper)
{
                       mylower=lower; 
                       myupper=upper;                
}


double UniformSampler::getnumber()
{
       int myrand=rand(); 
       while((myrand==0)||(myrand==RAND_MAX)){myrand = rand(); } //We want a number in (0, RAND_MAX). 
       
       double myuni = myrand/static_cast<double>(RAND_MAX); //Create a number in (0,1).
       return mylower+myuni*(myupper-mylower);        
}

//Exponential

ExponentialSampler::ExponentialSampler(double lambda)
{
 mylambda=lambda;                                              
} 

double ExponentialSampler::getnumber()
{
       UniformSampler myu; //since U(0, 1) is default option
       double tmp=-log(myu.getnumber());
       return tmp/mylambda;        
}

//Normal
NormalSampler::NormalSampler(double mu, double sigsqrd)
{
   mymu=mu; 
   mysigmasquared=sigsqrd;                                           
} 

double NormalSampler::getnumber()
{
       ExponentialSampler mye(0.5); 
       UniformSampler myu(0, 2*PI); 
       double standardnorm=sqrt(mye.getnumber())*cos(myu.getnumber()); 
       return mymu + sqrt(mysigmasquared)*standardnorm;       
}

//Doubleexp
DoubleExpSampler::DoubleExpSampler(void)
{
}

double DoubleExpSampler::getnumber()
{
  UniformSampler myu; //since U(0, 1) is default option
  double uninumber=myu.getnumber(); 
  double res; 
  if (uninumber > 0.5){res=-log(2.0*(1-uninumber));} 
  else {res=log(2.0*uninumber);}
  return res; 
}


//Functions for acceptance-rejection example
double stdnormaldens(double x)
{
       return exp(-0.5*x*x)/(sqrt(2.0*PI));       
}

double doubleexpdens(double x)
{
       return 0.5*exp(-fabs(x)); 
}

/////WARNING: CALL BY REFERENCE!
double specialRejectionSampling(int & numberrejections, double c)
{
       UniformSampler u; 
       DoubleExpSampler de; 
       bool accept=false; 
       double myx; 
       numberrejections=0; 
       while(accept==false){
                            myx = de.getnumber(); 
                            if (u.getnumber()*c*doubleexpdens(myx) <= stdnormaldens(myx)){accept=true; }
                            else {++numberrejections; } 
                                                 
       }
       return myx; 
       
       
}


