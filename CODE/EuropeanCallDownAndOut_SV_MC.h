#ifndef EUROPEANCALLDOWNANDOUT_SV_MC_H
#define EUROPEANCALLDOWNANDOUT_SV_MC_H

#include "BarrierPricer_SV.h"

class EuropeanCallDownAndOut_SV_MC : public BarrierPricer_SV
{
      public:
             EuropeanCallDownAndOut_SV_MC(int mc);
             virtual ~EuropeanCallDownAndOut_SV_MC(){};
             virtual double * getprice(double barr, double r, double alpha,
                 double beta, double gamma, double rho, double v0, double S0,
                 double K, double T, int m);
            
      private:
              int n;
};

#endif // EUROPEANCALLDOWNANDOUT_SV_MC_H
