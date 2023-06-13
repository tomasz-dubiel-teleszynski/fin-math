#ifndef BARRIERPRICER_SV_H
#define BARRIERPRICER_SV_H

class BarrierPricer_SV{
      
      public:
             virtual double * getprice(double barr, double r, double alpha,
                 double beta, double gamma, double rho, double v0, double S0,
                 double K, double T, int m) = 0;
};

#endif // BARRIERPRICER_SV_H


