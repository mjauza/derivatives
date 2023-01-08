#include "Chooser_option.h"
#include <cmath>
#include <algorithm>
#include "Call_option.h"

Chooser_option::Chooser_option(double K, double T) : c(K, T)
{
    this -> K = K;
    this -> T = T;
}

double Chooser_option::BS_Tc(double S_Tc, double Tc, double sigma, double r, double q)
{
    double call_bs = this -> c.BS(S_Tc,Tc,sigma,r,q);
    double v = call_bs + exp(-q*(T - Tc))* std::max(0.0, K*exp(-(r-q)*(T-Tc)) - S_Tc);
    return v;
}

Chooser_option::~Chooser_option()
{
    //dtor
}
