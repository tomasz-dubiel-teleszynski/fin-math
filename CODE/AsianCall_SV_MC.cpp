#include <algorithm> // max
#include <cmath> // exp, sqrt, log
#include <new> // new

#include "AsianCall_SV_MC.h"
#include "Sampler.h" // normal sampler

using namespace std;

AsianCall_SV_MC::AsianCall_SV_MC(int mc)
{
       n = mc;
}

double * AsianCall_SV_MC::getprice(double r, double alpha, double beta,
    double gamma, double rho, double v0, double S0, double K, double T, int m)
{      
       double eps = 1e-6;
       double h = T/m;
          
       NormalSampler NS;
       double Zv;
       double Zs;
       
       double *v = new double[m+1];
       double *S = new double[m+1];
       double *f = new double[n];
       
       v[0] = v0;
       S[0] = log(S0);
       double F = 0;
       
       // ---------
       double sumS;
       // ---------
       
       for(int i = 0; i < n; i++)
       {  
               
               // ------
               sumS = 0;
               // ------               
               
               for(int j = 1; j <=m; j++)
               {
                       Zv = NS.getnumber();
                       Zs = rho*Zv + sqrt(1-pow(rho,2))*NS.getnumber();
                       
                       v[j] = v[j-1] + alpha * ( beta - max(eps,v[j-1]) ) * h +
                           gamma * sqrt( max(eps,v[j-1]) ) * Zv * sqrt(h);
                       S[j] = S[j-1] + ( r - 0.5 * max(eps,v[j-1]) ) * h +
                           sqrt( max(eps,v[j-1]) ) * Zs * sqrt(h);
                       
                       // ----------
                       sumS += exp(S[j]);
                       // ----------                                        
               }
               f[i] = exp(-r*T) * max( sumS/m - K, 0.0);
               F += f[i];       
       }
       F = F/n;
       
       double sn = 0;
       for(int i = 0; i < n; i++)
       {
               sn += (f[i]-F)*(f[i]-F);       
       }
       sn = sn/n;
       
       double a = 1.96;
       double CIl = F-a*sn/sqrt(n);
       double CIu = F+a*sn/sqrt(n);
       
       static double result[2];
       result[0] = CIl;
       result[1] = F;
       result[2] = CIu;
       
       delete [] v;
       delete [] S;
       delete [] f;
       
       return result;
}
