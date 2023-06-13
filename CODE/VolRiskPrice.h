#ifndef VOLRISKPRICE_H
#define VOLRISKPRICE_H

#include "VolRiskClass.h"

class VolRiskPrice : public VolRiskClass
{
      public:
            VolRiskPrice(double _r, double _alpha, double _beta, double _gamma,
                double _rho,double _v0,double _S0,double _K,double _T);
            virtual ~VolRiskPrice(){};
            virtual double Value(double lambda);
      
      private:
             double S0;
             double K;
             double T;
             double r;
             
             double v0;
             double alpha;
             double beta;
             double gamma;
             double rho;
};

#endif // VOLRISKPRICE_H
