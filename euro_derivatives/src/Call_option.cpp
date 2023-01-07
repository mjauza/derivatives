#include "Call_option.h"
#include <cmath>
Call_option::Call_option(double K, double T)
{
    this -> K = K;
    this -> T = T;
}

double Call_option::normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x/std::sqrt(2))/2;
}

double Call_option::BS(double S, double t, double sigma, double r)
{
    double d1 = 1/(sigma * sqrt(T-t)) * ( log(S/K) + (r+ pow(sigma,2)/2)*(T-t) );
    double d2 = d1 - sigma*sqrt(T-t);
    double C = normalCDF(d1)*S - normalCDF(d2)*K*exp(-r*(T - t));
    return C;
}

Call_option::~Call_option()
{
    //dtor
}
