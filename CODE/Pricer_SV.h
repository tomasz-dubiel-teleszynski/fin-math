#ifndef PRICER_SV_H
#define PRICER_SV_H

class Pricer_SV{
      
      public:
             virtual double * getprice(double r, double alpha, double beta,
                 double gamma, double rho, double v0, double S0, double K,
                 double T, int m) = 0;

};

#endif // PRICER_SV_H


