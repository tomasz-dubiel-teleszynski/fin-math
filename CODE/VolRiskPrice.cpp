#include "VolRiskPrice.h"

// declare function returning closed form price from Heston model
double getpriceCF(double r, double alpha, double beta, double gamma,
    double rho, double lambda, double v0, double S0, double K, double T);

using namespace std;

// constructor
VolRiskPrice::VolRiskPrice(double _r, double _alpha, double _beta, double _gamma,
    double _rho,double _v0,double _S0,double _K,double _T){
    
    S0 = _S0;
    K  = _K;
    T  = _T;
    r  = _r;
    
    v0    = _v0;
    alpha = _alpha;
    beta  = _beta;
    gamma = _gamma;
    rho   = _rho;
    
}

// define member function
double VolRiskPrice::Value(double lambda)
{
    return getpriceCF(r,alpha,beta,gamma,rho,lambda,v0,S0,K,T);       
}
