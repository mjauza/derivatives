#include "Put_option.h"
#include <cmath>
#define _USE_MATH_DEFINES
Put_option::Put_option(double K, double T)
{
    this -> K = K;
    this -> T = T;
}
double Put_option::normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x/std::sqrt(2))/2;
}

double Put_option::BS(double S, double t, double sigma, double r, double q)
{
    double d1 = 1/(sigma * sqrt(T-t)) * ( log(S/K) + (r - q + pow(sigma,2)/2)*(T-t) );
    double d2 = d1 - sigma*sqrt(T-t);
    double P = normalCDF(-d2)*K*exp(-r*(T - t)) - normalCDF(-d1)*S*exp(-q*(T-t));
    return P;
}

double Put_option::delta(double S, double t, double sigma, double r)
{
    double d1 = 1/(sigma * sqrt(T-t)) * ( log(S/K) + (r+ pow(sigma,2)/2)*(T-t) );
    return -normalCDF(-d1);
}

double Put_option::theta(double S, double t, double sigma, double r)
{
    double d1 = 1/(sigma * sqrt(T-t)) * ( log(S/K) + (r+ pow(sigma,2)/2)*(T-t) );
    double d2 = d1 - sigma*sqrt(T-t);
    double tht = -(1/sqrt(2*M_PI))*(S * exp(-pow(d1,2)/2) * sigma / (2*sqrt(T-t))) + r*K*exp(-r*(T-t))*normalCDF(-d2);
    return tht;
}

double Put_option::gamma(double S, double t, double sigma, double r)
{
    double d1 = 1/(sigma * sqrt(T-t)) * ( log(S/K) + (r+ pow(sigma,2)/2)*(T-t) );
    double gam = exp(-pow(d1,2))/(S*sigma*sqrt(2*M_PI*(T-t)));
    return gam;
}

double Put_option::vega(double S, double t, double sigma, double r)
{
    double d1 = 1/(sigma * sqrt(T-t)) * ( log(S/K) + (r+ pow(sigma,2)/2)*(T-t) );
    double v = S*sqrt(T-t)*exp(-pow(d1,2)/2)/sqrt(2*M_PI);
    return v;
}

double Put_option::rho(double S, double t, double sigma, double r)
{
    double d1 = 1/(sigma * sqrt(T-t)) * ( log(S/K) + (r+ pow(sigma,2)/2)*(T-t) );
    double d2 = d1 - sigma*sqrt(T-t);
    double rh = -(T-t)*K*exp(-r*(T-t))*normalCDF(-d2);
    return rh;
}

Put_option::~Put_option()
{
    //dtor
}