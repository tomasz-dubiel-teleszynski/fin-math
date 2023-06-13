////////////////////////////////////////////////////

// paths at the LSE --------------------------------

//Q:/A20010_BLOODSHED_4992/Dev-Cpp/lib/libgsl.a
//Q:/A20010_BLOODSHED_4992/Dev-Cpp/lib/libgslcblas.a

////////////////////////////////////////////////////

////////////////////////////////////////////////////

// paths  at home ----------------------------------

//C:/Dev-Cpp/lib/libgsl.a
//C:/Dev-Cpp/lib/libgslcblas.a

////////////////////////////////////////////////////

#include <cstdlib>              // EXIT_SUCCESS, EXIT_FAILURE
#include <vector>               // vector
#include <iostream>             // cout
#include <stdio.h>              // writing output file
#include <cmath>                // exp, fabs
#include <gsl/gsl_multimin.h>   // Simplex algorithm of Nelder and Mead
#include "EuropeanCall_SV_MC.h" // European call option under Heston stochastic 
                                // volatility model priced using Monte Carlo
#include "EuropeanCallDownAndOut_SV_MC.h" // Down-and-out European call option 
                                // under Heston stochastic volatility model 
                                // priced using Monte Carlo 
#include "AsianCall_SV_MC.h"    // Asian call option under Heston stochastic 
                                // volatility model priced using Monte Carlo 
#include "VolRiskClass.h"       // volatility of risk base class
#include "VolRiskPrice.h"       // volatility of risk price derived class

using namespace std;

// declare functions used here which are defined in functions.cpp
int minf(vector<double>& calvals, double *fpar, double *stvals);
int Feller(double alpha, double beta, double gamma);
double rmsef(double *prices, double *pricesCAL, int n);
double Bisection(VolRiskClass *FCptr, double a, double b, double c, double eps);

////////////////////////////////////////////////////////////////////////////////
// TO TEST ENTIRE CODE DECREASE "m" (FIRST VARIABLE BELOW) TO "m = 100" FOR ////
// SPEED, HOWEVER IT WILL NOT AFFECT SPEED OF NUMERICAL INTEGRATION WITHIN /////
// MINIMIZATION ONLY SIMULATING PATHS THUS OPTIMIZATION TIME REMAINS SAME //////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// MAIN PROGRAM ////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main(void)
{
  //////////////////////////////////////////////////////////////////////////////
  // PROBLEM 1 /////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  // declare and specify Euler scheme and Monte Carlo simulation parameters
  // number of Euler scheme steps
  int m = 1000;   
  // number of Monte Carlo simulations
  int n = 100000;
  
  // declare and specify general parameters for pricing European call option
  double r  =  0.0;
  double S0 = 10.0;
  double K  = 10.0;
  double T  =  1.0;  // equivalent to 1 year
  
  // declare and specify parameters for Heston stochastic volatility model
  double v0    =  0.04;
  double alpha =  2.0; 
  double beta  =  0.04;
  double gamma =  0.2; 
  double rho   = -0.2;
  // check Feller condition for above start values        
  if(Feller(alpha,beta,gamma))
  {
      cout << "Feller condition not met by specified parameters!" <<  endl;
      system("PAUSE");
      return EXIT_FAILURE;
  }
  
  // instantiate European call option under stochastic volatility using Monte 
  // Carlo pricer class
  EuropeanCall_SV_MC ecmcsv(n);
  
  // European call option price
  double *ecmcsvp = ecmcsv.getprice(r,alpha,beta,gamma,rho,v0,S0,K,T,m); 
  
  // open file
  FILE *pricingEC = fopen("PricingEC.txt","w");
  // screen output
  printf("European call option price:\n");
  printf("LCIB: %.4f PRICE: %.4f UCIB: %.4f\n",ecmcsvp[0],ecmcsvp[1],
      ecmcsvp[2]);
  // file output
  fprintf(pricingEC,"European call option price:\n");
  fprintf(pricingEC,"LCIB: %.4f PRICE: %.4f UCIB: %.4f\n",ecmcsvp[0],ecmcsvp[1],
      ecmcsvp[2]);
  // close file
  fclose(pricingEC);
  
  //////////////////////////////////////////////////////////////////////////////
  // PROBLEM 3 /////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  // problem resultant from the fact that later we used closed form solution
  // to Heston stochastic volatility model for European call option price
  // hence we need to calculate price of volatility risk (volatility risk
  // premium) using bisection method and reference price obtained in problem 1
  double refprice = ecmcsvp[1];
  // VolRiskPrice class instantiation
  VolRiskPrice vrp(r,alpha,beta,gamma,rho,v0,S0,K,T);
  VolRiskPrice *VRPptr = &vrp;
  // specify tolerance
  double eps = 0.0000000000000001;
  // specify bounds for bisection method
  double lbound = -5;
  double ubound =  5;
  // calculate volatility risk premium using bisection method
  double lambda =  Bisection(VRPptr,lbound,ubound,refprice,eps);   
  
  // open file 
  FILE *volriskpr = fopen("VolRiskPremium.txt","w");
  // terminal output ----------------------
  printf("Price of volatility risk (volatility risk premium): %.4f\n",lambda);
  // file output
  fprintf(volriskpr,"Price of volatility risk (volatility risk premium): %.4f\n"
      ,lambda);   
  // close file
  fclose(volriskpr);
  
  // specify initial stock value for pricing European call option which is 
  // fixed in minimization and differs from this specified above in problem 1
  S0 = 14.97; 
  
  // specify minimization criterion (1 for SSE and 0 for SAD)
  double cr = 1;
  
  // collect general parameters for pricing European call option which are 
  // fixed in minimization and minimization criterion in double array
  double fpar[7] = {T, S0, r, (double)m, (double)n, lambda, cr};  
  
  // declare and specify start values for Heston stochastic volatility model,
  // apply logarithmic transormation due to exponentiation in objective function
  // for convenience choose corresponding values from problem 1
  double v00    =  log(v0);
  double alpha0 =  log(alpha);
  double beta0  =  log(beta);
  double gamma0 =  log(gamma);
  double rho0   = -log(2.0/(rho+1.0)-1);
  // check Feller condition for above start values
  if(Feller(exp(alpha0), exp(beta0), exp(gamma0)))
  {
      cout << "Feller condition not met by start values!" <<  endl;
      system("PAUSE");
      return EXIT_FAILURE;
  }
  
  //collect above start values in a double array
  double stvals[5] = {v00,alpha0,beta0,gamma0,rho0};
  
  // calibrate to MARKET DATA using SSE criterion
  // define double vector to store SSE calibration results
  vector<double> calvalsSSE(5);
  // begin SSE minimization
  int retSSE = minf(calvalsSSE,fpar,stvals);
  // check status of SSE minimization
  if(retSSE)
  {
     cout << "SSE minimization failed!" <<  endl;
     system("PAUSE");
     return EXIT_FAILURE;
  } 
  
  // calibrate to MARKET DATA using SAD criterion
  // change minimization criterion (1 for SSE and 0 for SAD)
  cr = 0; 
  // inset change from SSE to SAD into double array of fixed parameters
  fpar[6] = cr;
  // define double vector to store SAD calibration results
  vector<double> calvalsSAD(5);
  // begin SAD minimization
  int retSAD = minf(calvalsSAD,fpar,stvals);
  // check status of SAD minimization
  if(retSAD)
  {
     cout << "SAD minimization failed!" <<  endl;
     system("PAUSE");
     return EXIT_FAILURE;
  }
  
  // calculate post-SSE-calibration European call option prices    
  // assign SSE-calibrated parameter values to variables
  double _v0    = calvalsSSE[0];
  double _alpha = calvalsSSE[1];
  double _beta  = calvalsSSE[2];
  double _gamma = calvalsSSE[3];
  double _rho   = calvalsSSE[4];
  // check Feller condition for above SSE-calibrated parameter values
  if(Feller(_alpha, _beta, _gamma))
  {
      cout << "Feller condition not met by SSE-calibrated values!" <<  endl;
      system("PAUSE");
      return EXIT_FAILURE;
  }
  
  // MARKET DATA ----------------------------------------------------
  // constant number of market data observations
  const int nStrikes = 7;
  // given strikes
  double strikes[nStrikes] = {10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0};
  // given market prices
  double prices[nStrikes] = {5.43, 3.35, 1.85, 0.9, 0.4, 0.18, 0.09};
  // ----------------------------------------------------------------
  
  // define double array to store post-SSE-calibration European call
  // option prices
  double pricesSSE[nStrikes];
  
  // open file 
  FILE *pricingECSSE = fopen("PricingECSSE.txt","w");
  // terminal output ----------------------
  printf("Post-SSE-calibration European call option prices:\n");
  // file output
  fprintf(pricingECSSE,"Post-SSE-calibration European call option prices:\n");   
  // loop over strikes
  for(int i = 0; i<nStrikes; i++)
  {
      // post-SSE-calibration European call option price
      ecmcsvp = ecmcsv.getprice(r,_alpha,_beta,_gamma,_rho,_v0,S0,
          strikes[i],T,m);
      pricesSSE[i] = ecmcsvp[1];
            
      // terminal output
      printf("LCIB: %.4f PRICE: %.4f UCIB: %.4f\n",ecmcsvp[0],ecmcsvp[1],
          ecmcsvp[2]);
      // file output
      fprintf(pricingECSSE,"LCIB: %.4f PRICE: %.4f UCIB: %.4f\n",ecmcsvp[0],
          ecmcsvp[1],ecmcsvp[2]);
  }
  // calculate RMSE for post-SSE-calibration prices
  double rmseSSE = rmsef(prices, pricesSSE, nStrikes);
   // terminal output
  printf(
      "RMSE for market and post-SSE-calibration European call option prices: %.4f\n",
      rmseSSE);
  // file output
  fprintf(pricingECSSE,
      "RMSE for market and post-SSE-calibration European call option prices: %.4f\n",
      rmseSSE);   
  // close file
  fclose(pricingECSSE);
  
  // calculate post-SAD-calibration European call option prices 
  // assign SAD-calibrated parameter values to variables
  double _v0_    = calvalsSAD[0];
  double _alpha_ = calvalsSAD[1];
  double _beta_  = calvalsSAD[2];
  double _gamma_ = calvalsSAD[3];
  double _rho_   = calvalsSAD[4];
  // check Feller condition for above SAD-calibrated parameter values
  if(Feller(_alpha_,_beta_,_gamma_))
  {
      cout << "Feller condition not met by SAD-calibrated values!" <<  endl;
      system("PAUSE");
      return EXIT_FAILURE;
  }
  
  // define double array to store post-SAD-calibration European call
  // option prices
  double pricesSAD[nStrikes];
  
  // open file 
  FILE *pricingECSAD = fopen("PricingECSAD.txt","w");
  // terminal output
  printf("Post-SAD-calibration European call option prices:\n");
  // file output
  fprintf(pricingECSAD,"Post-SAD-calibration European call option prices:\n");  
  // loop over strikes
  for(int i = 0; i<nStrikes; i++)
  {
      // post-SAD-calibration European call option price
      ecmcsvp = ecmcsv.getprice(r,_alpha_,_beta_,_gamma_,_rho_,_v0_,S0,
          strikes[i],T,m);
      pricesSAD[i] = ecmcsvp[1];
      
      // terminal output
      printf("LCIB: %.4f PRICE: %.4f UCIB: %.4f\n",ecmcsvp[0],ecmcsvp[1],
          ecmcsvp[2]);
      // file output
      fprintf(pricingECSAD,"LCIB: %.4f PRICE: %.4f UCIB: %.4f\n",ecmcsvp[0],
          ecmcsvp[1],ecmcsvp[2]);
  }
   // calculate RMSE for post-SAD-calibration prices
  double rmseSAD = rmsef(prices, pricesSAD, nStrikes);
   // terminal output
  printf(
      "RMSE for market and post-SAD-calibration European call option prices: %.4f\n",
      rmseSAD);
  // file output
  fprintf(pricingECSAD,
      "RMSE for market and post-SAD-calibration European call option prices: %.4f\n",
      rmseSAD);   
  // close file
  fclose(pricingECSAD);
  
  //////////////////////////////////////////////////////////////////////////////
  // PROBLEM 4 /////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  // calculate post-SSE-calibration down-and-out European call option prices    
  // instantiate down-and-out European call option under stochastic volatility 
  // using Monte Carlo pricer class
  EuropeanCallDownAndOut_SV_MC ecdaomcsv(n);
  // define pointer to price from above pricer
  double *ecdaomcsvp;
  // define price barrier
  double barr = 9;   
  
  // define double array to store post-SSE-calibration down-and-out European
  // call option prices
  double pricesSSEDAO[nStrikes];
    
  // open file 
  FILE *pricingSSEDAO = fopen("PricingSSEDAO.txt","w");
  // terminal output
  printf("Post-SSE-calibration down-and-out European call option prices:\n");
  // file output
  fprintf(pricingSSEDAO,
      "Post-SSE-calibration down-and-out European call option prices:\n");  
  // loop over strikes
  for(int i = 0; i<nStrikes; i++)
  {
      // post-SSE-calibration down-and-out European call option price
      ecdaomcsvp = ecdaomcsv.getprice(barr,r,_alpha,_beta,_gamma,_rho,_v0,S0,
          strikes[i],T,m);
      pricesSSEDAO[i] = ecdaomcsvp[1];
      
      // terminal output
      printf("LCIB: %.4f PRICE: %.4f UCIB: %.4f\n",ecdaomcsvp[0],ecdaomcsvp[1],
          ecdaomcsvp[2]);
      // file output
      fprintf(pricingSSEDAO,"LCIB: %.4f PRICE: %.4f UCIB: %.4f\n",ecdaomcsvp[0],
          ecdaomcsvp[1],ecdaomcsvp[2]);
  } 
  // close file
  fclose(pricingSSEDAO);
  
  // calculate post-SAD-calibration down-and-out European call option prices  
  // define double array to store post-SAD-calibration down-and-out European
  // call option prices
  double pricesSADDAO[nStrikes];
  
   // open file 
  FILE *pricingSADDAO = fopen("PricingSADDAO.txt","w");
  // terminal output
  printf("Post-SAD-calibration down-and-out European call option prices:\n");
  // file output
  fprintf(pricingSADDAO,
      "Post-SAD-calibration down-and-out European call option prices:\n");  
  // loop over strikes
  for(int i = 0; i<nStrikes; i++)
  {
      // post-SAD-calibration down-and-out European call option price
      ecdaomcsvp = ecdaomcsv.getprice(barr,r,_alpha,_beta,_gamma,_rho,_v0,S0,
          strikes[i],T,m);
      pricesSADDAO[i] = ecdaomcsvp[1];
      
      // terminal output
      printf("LCIB: %.4f PRICE: %.4f UCIB: %.4f\n",ecdaomcsvp[0],ecdaomcsvp[1],
          ecdaomcsvp[2]);
      // file output
      fprintf(pricingSADDAO,"LCIB: %.4f PRICE: %.4f UCIB: %.4f\n",
          ecdaomcsvp[0],ecdaomcsvp[1],ecdaomcsvp[2]);
  } 
  // close file
  fclose(pricingSADDAO);
  
  // calculate post-SSE-calibration Asian call option prices    
  // instantiate Asian call option under stochastic volatility using Monte 
  // Carlo pricer class
  AsianCall_SV_MC acmcsv(n);
  // define pointer to price from above pricer
  double *acmcsvp;
  
  // define double array to store post-SSE-calibration Asian call option prices
  double pricesSSEAC[nStrikes];
  
  // open file 
  FILE *pricingSSEAC = fopen("PricingSSEAC.txt","w");
  // terminal output
  printf("Post-SSE-calibration Asian call option prices:\n");
  // file output
  fprintf(pricingSSEAC,"Post-SSE-calibration Asian call option prices:\n");  
  // loop over strikes
  for(int i = 0; i<nStrikes; i++)
  {
      // post-SSE-calibration Asian call option price
      acmcsvp = acmcsv.getprice(r,_alpha,_beta,_gamma,_rho,_v0,S0,
          strikes[i],T,m);
      pricesSSEAC[i] = acmcsvp[1];
      
      // terminal output
      printf("LCIB: %.4f PRICE: %.4f UCIB: %.4f\n",acmcsvp[0],acmcsvp[1],
          acmcsvp[2]);
      // file output
      fprintf(pricingSSEAC,"LCIB: %.4f PRICE: %.4f UCIB: %.4f\n",acmcsvp[0],
          acmcsvp[1],acmcsvp[2]);
  }
  // close file
  fclose(pricingSSEAC);
  
  // calculate post-SAD-calibration Asian call option prices   
  // define double array to store post-SAD-calibration Asian call option prices
  double pricesSADAC[nStrikes];
  
  // open file 
  FILE *pricingSADAC = fopen("PricingSADAC.txt","w");
  // terminal output
  printf("Post-SAD-calibration Asian call option prices:\n");
  // file output
  fprintf(pricingSADAC,"Post-SAD-calibration Asian call option prices:\n");  
  // loop over strikes
  for(int i = 0; i<nStrikes; i++)
  {
      // post-SAD-calibration Asian call option price
      acmcsvp = acmcsv.getprice(r,_alpha_,_beta_,_gamma_,_rho_,_v0_,S0,
          strikes[i],T,m);
      pricesSADAC[i] = acmcsvp[1];
      
      // terminal output
      printf("LCIB: %.4f PRICE: %.4f UCIB: %.4f\n",acmcsvp[0],acmcsvp[1],
          acmcsvp[2]);
      // file output
      fprintf(pricingSADAC,"LCIB: %.4f PRICE: %.4f UCIB: %.4f\n",acmcsvp[0],
          acmcsvp[1],acmcsvp[2]);
  }
  // close file
  fclose(pricingSADAC);
  
  system("PAUSE");
  return EXIT_SUCCESS;
}
