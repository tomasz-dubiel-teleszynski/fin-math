#ifndef ASIANCALL_SV_MC_H
#define ASIANCALL_SV_MC_H

#include "Pricer_SV.h"

class AsianCall_SV_MC : public Pricer_SV
{
      public:
             AsianCall_SV_MC(int mc);
             virtual ~AsianCall_SV_MC(){};
             virtual double * getprice(double r, double alpha, double beta,
                 double gamma, double rho, double v0, double S0, double K,
                 double T, int m);
            
      private:
              int n;
};

#endif // ASIANCALL_SV_MC_H
