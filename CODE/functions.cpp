#include <complex> // complex
#include <vector> // vector
#include <cmath> // exp, fabs
#include <gsl/gsl_multimin.h> // Simplex algorithm of Nelder and Mead
#include "EuropeanCall_SV_MC.h" // European call option under stochastic 
                                // volatility priced using Monte Carlo 
                                // simulation
#include "VolRiskPrice.h"       // volatility of risk price class

using namespace std;

// declare functions called here by functions used then in main.cpp
double realf(double phi, double alpha, double beta, double gamma,
    double rho, double lambda, double v0, double S0, double K, double r, double T,
    int fNum);
double trapz(vector<double> X, vector<double> Y);
double objf (const gsl_vector *v, void *params);
double getpriceCF(double r, double alpha, double beta, double gamma,
    double rho, double lambda, double v0, double S0, double K, double T);
double Bisection(VolRiskClass *FCptr, double a, double b, double c, double epsilon);

// -----------------------------------------------------------------------------
// Bisection method to calculate price of volatility risk ----------------------
// -----------------------------------------------------------------------------
double Bisection(VolRiskClass *FCptr, double a, double b, double c, double eps)
{
    double l = a;
    double r = b;
    while((r-l)>eps)
    {
        if( ( FCptr->Value(l)-c ) * ( FCptr->Value(0.5*(l+r))-c ) <= 0 )
        {
            l = l;
            r = 0.5*(l+r);
        }
        else
        {
            l = 0.5*(l+r);
            r = r;
        }
    }
    return l;
}

// -----------------------------------------------------------------------------
// European call option price from Heston model in closed form------------------
// -----------------------------------------------------------------------------
double getpriceCF(double r, double alpha, double beta, double gamma,
    double rho, double lambda, double v0, double S0, double K, double T)
{
    const double pi = 3.141592653589793238462643;	
	
    // create partition	
	double infty = 100.0;
	double step  =   0.01;
    int nPoints  = int(infty/step);
	
	// declare complex vectors
	vector<double> PHI(nPoints);
	vector<double>  F1(nPoints);
	vector<double>  F2(nPoints);
	
	// grid evaluation
    for (int j=0; j<=nPoints-1; j++)
    {
		if(j==0)
        {
            PHI[j] = 0.0000001;
        }
        else
        {
		    PHI[j]  = double(j)*step;
        }
        F1[j] = realf(PHI[j], alpha, beta, gamma, rho, lambda, v0, S0, K, r, T, 1);
		F2[j] = realf(PHI[j], alpha, beta, gamma, rho, lambda, v0, S0, K, r, T, 2);
	}
	
	// integration
    double int1 = trapz(PHI, F1);
	double int2 = trapz(PHI, F2);
	
	// probabilities
    double P1 = 0.5 + int1/pi;
	double P2 = 0.5 + int2/pi;
	
	// just in case
    if(P1<0){P1=0;}
	if(P1>1){P1=1;}
    if(P2<0){P2=0;}
    if(P2>1){P2=1;}
	
	// return Heston price in closed form
    return S0*P1 - K*exp(-r*T)*P2;
}

// -----------------------------------------------------------------------------
// real part of characteristic function ----------------------------------------
// -----------------------------------------------------------------------------
double realf(double phi, double alpha, double beta, double gamma,
    double rho, double lambda, double v0, double S0, double K, double r, double T,
    int fNum)
{
    // create complex 1 and i
    complex<double>  one(1.0, 0.0);
	complex<double>    i(0.0, 1.0);
	
	// declare double variables
    double a = alpha*beta;
    double u;
    double b;
        
    // declare complex variables
	complex<double> d;
	complex<double> g;
	complex<double> C;
	complex<double> D;
	
	// make distinction
	if (fNum==1)
	{
        u = 0.5;
        b = alpha + lambda - rho*gamma;           
    }
	else
	{
	    u = -0.5;    
        b = alpha + lambda;
    }

    // calculate elements of characteristic function
    d = sqrt( pow(rho*gamma*phi*i-b,2) - pow(gamma,2)*( 2*u*phi*i - pow(phi,2) ) ); 
       
    g = (b-rho*gamma*phi*i+d)/(b-rho*gamma*phi*i-d);
     
    C = r*phi*i*T + alpha*beta/pow(gamma,2)*( (b-rho*gamma*phi*i+d)*T - 
        log( ((one-g*exp(d*T))/(one-g))*((one-g*exp(d*T))/(one-g)) ) );
        
    D = (b-rho*gamma*phi*i+d)/pow(gamma,2)*((one-exp(d*T))/(one-g*exp(d*T)));
    
    // return real part
    return real( exp( -i*phi*log(K) ) * exp( C + D*v0 + i*phi*log(S0) ) / (i*phi) );
}

// -----------------------------------------------------------------------------
// numerical integration by trapeizoidal rule ----------------------------------
// -----------------------------------------------------------------------------
double trapz(vector<double> X, vector<double> Y)
{
    int n = X.size();
	double sum = 0;
	for (int i=1; i<=n-1; i++)
		sum += 0.5*(X[i] - X[i-1])*(Y[i-1] + Y[i]);
	return sum;
}

// -----------------------------------------------------------------------------
// root mean-squared error -----------------------------------------------------
// -----------------------------------------------------------------------------
double rmsef(double * prices, double * pricesCAL, int n)
{
    double rmse = 0.0;
    for(int i = 0; i<n; i++)
    {
        rmse += (prices[i]-pricesCAL[i])*(prices[i]-pricesCAL[i]);
    }
    rmse = rmse/n;
    return sqrt(rmse);
}

// -----------------------------------------------------------------------------
// Feller condition check ------------------------------------------------------
// -----------------------------------------------------------------------------
int Feller(double alpha, double beta, double gamma)
{
    if(2*alpha*beta-gamma*gamma<0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

// -----------------------------------------------------------------------------
// minimization objective function ---------------------------------------------
// -----------------------------------------------------------------------------
double objf(const gsl_vector *v, void *params)
{  
  // MARKET DATA ---------------------------------------------
  // given strikes ///////////////////////////////////////////
  double strikes[7] = {10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0};
  // given market prices /////////////////////////////////////
  double prices[7] = {5.43, 3.35, 1.85, 0.9, 0.4, 0.18, 0.09};
  
  // cast void to double
  double *par = (double *)params;
  
  // allocate fixed parameter values to variables
  double T      =      par[0];
  double S0     =      par[1];
  double r      =      par[2];
     int m      = (int)par[3];
     int n      = (int)par[4];
  double lambda =      par[5];   
     int cr     = (int)par[6];
  
  // allocate iterated parameter valies to variables
  double v0, alpha, beta, gamma, rho;
  v0    = exp(gsl_vector_get(v, 0));
  alpha = exp(gsl_vector_get(v, 1));
  beta  = exp(gsl_vector_get(v, 2));
  gamma = exp(gsl_vector_get(v, 3));
  rho   = 2.0/(1.0+exp(-gsl_vector_get(v, 4)))-1.0;
  
  double p;
  double crit = 0;
  for(int i =0; i<7; i++)
  {
      // calculate  closed form European call option price
      p = getpriceCF(r,alpha,beta,gamma,rho,lambda,v0,S0,strikes[i],T);
      
      if(cr)
      {
          // use SSE criterion
          crit += (prices[i]-p)*(prices[i]-p);
      }
      else
      {
          // use SAD criterion
          crit += fabs(prices[i]-p);
      }
  }

  // penalty term multiplier
  double penmult = 1e6;
  
  // penalty term
  double penalty = max(0.0,gamma*gamma-2.0*alpha*beta);
  
  return crit + penmult * penalty; 
}

// -----------------------------------------------------------------------------
// minimizing procedure --------------------------------------------------------
// -----------------------------------------------------------------------------
int minf(vector<double>& calvals, double * fpar, double * stvals)
{        
    // Simplex algorithm of Nelder and Mead instantiation
    const gsl_multimin_fminimizer_type *SNM = gsl_multimin_fminimizer_nmsimplex;     
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;
    
    // allocate start values to variables
    double v00    = stvals[0];
    double alpha0 = stvals[1];
    double beta0  = stvals[2];
    double gamma0 = stvals[3];
    double rho0   = stvals[4];
    
    // collect variables with start values into GSL vector                
    x = gsl_vector_alloc(5);
    gsl_vector_set(x,0,v00); 
    gsl_vector_set(x,1,alpha0);
    gsl_vector_set(x,2,beta0);  
    gsl_vector_set(x,3,gamma0);
    gsl_vector_set(x,4,rho0);  

    // set initial step size to 1 for all
    ss = gsl_vector_alloc(5);
    gsl_vector_set_all(ss,1);
  
    // set up minimizer 
    minex_func.n = 5;
    minex_func.f = objf;
    minex_func.params = fpar;
    s = gsl_multimin_fminimizer_alloc(SNM, 5);
    gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

    // create output file ----------------------------
    FILE *calibration;
    if(fpar[6])
    {
        calibration = fopen("CalibrationSSE.txt","w");
        // terminal output --------------
        printf("SSE calibration:\n");
        // file output
        fprintf(calibration,"SSE calibration:\n");
        // ------------------------------
    }
    else
    {
        calibration = fopen("CalibrationSAD.txt","w");
        // terminal output --------------
        printf("SAD calibration:\n");
        // file output
        fprintf(calibration,"SAD calibration:\n");
        // ------------------------------
    }
    // -----------------------------------------------
  
    // auxiliaries
    size_t iter = 0;
    int status;
    double size;
    
    do
    {
        // iterating
        iter++;
      
        // check status
        status = gsl_multimin_fminimizer_iterate(s);
        if (status) 
            break;

        // check termination criterion
        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, 1e-6);
        
        // check status again
        if (status == GSL_SUCCESS)
        {
            // terminal output -------------------------
            printf ("converged to minimum at\n");
            // file output
            fprintf(calibration,"converged to minimum at\n");
            // -----------------------------------------
        }

        // terminal output ---------------------------------------------
        printf ("%5d %.4f %.4f %.4f %.4f %.4f f() = %.6f size = %.6f\n",
               iter,
               exp(gsl_vector_get(s->x, 0)), 
               exp(gsl_vector_get(s->x, 1)),
               exp(gsl_vector_get(s->x, 2)),
               exp(gsl_vector_get(s->x, 3)),
               2.0/(1+exp(-gsl_vector_get(s->x, 4)))-1,
               s->fval,
               size);
        // file output
        fprintf(calibration,
                "%5d %.4f %.4f %.4f %.4f %.4f f() = %.6f size = %.6f\n", 
                iter,
                exp(gsl_vector_get(s->x, 0)), 
                exp(gsl_vector_get(s->x, 1)),
                exp(gsl_vector_get(s->x, 2)),
                exp(gsl_vector_get(s->x, 3)),
                2.0/(1+exp(-gsl_vector_get(s->x, 4)))-1,
                s->fval,
                size);
        // ------------------------------------------------------------
    } // set maximum number of iterations below to 2000   
    while (status == GSL_CONTINUE && iter < 2001);
   
    // close output file -
    fclose(calibration);
    // -------------------
    
    // assign calibrated parameter values to output vector
    calvals.at(0) = exp(gsl_vector_get(s->x, 0));             // v0
    calvals.at(1) = exp(gsl_vector_get(s->x, 1));             // alpha
    calvals.at(2) = exp(gsl_vector_get(s->x, 2));             // beta
    calvals.at(3) = exp(gsl_vector_get(s->x, 3));             // gamma;
    calvals.at(4) = 2.0/(1+exp(-gsl_vector_get(s->x, 4)))-1;  // rho
    
    // release GSL vectors and minimizer
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
    
    // report final status
    return status;
}
