#include "Call_option.h"
#include <cmath>
#define _USE_MATH_DEFINES
#include <iostream>
#include "Simulation_methods.h"
#include <numeric>
Call_option::Call_option(double K, double T)
{
    this -> K = K;
    this -> T = T;
}

double Call_option::normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x/std::sqrt(2))/2;
}

double Call_option::BS(double S, double t, double sigma, double r, double q)
{
    double d1 = 1/(sigma * sqrt(T-t)) * ( log(S/K) + (r - q + pow(sigma,2)/2)*(T-t) );
    double d2 = d1 - sigma*sqrt(T-t);
    double C = normalCDF(d1)*S*exp(-q*(T-t)) - normalCDF(d2)*K*exp(-r*(T - t));
    return C;
}

double Call_option::MC_BS(double S, double t, double sigma, double r, int N){

    // simulate ST
    std::vector<double> ST = Simulation_methods::geom_BM_t(N, T-t, r, sigma, S);

    // calculate CT
    std::vector<double> CT;
    for(int i = 0; i<N; i++){
        double c = std::max(ST[i] - K, 0.0);
        CT.push_back(c);
    }

    /*
    for(double s : ST){
        double c = std::max(s - K, 0.0);
        std::cout << c << "\n";
    }
    */



    // avergae
    double C_sum = std::accumulate(CT.begin(), CT.end(), 0);
    double C_avg = C_sum / N;

    // discount
    double Ct = exp(-r*(T-t))*C_avg;
    return Ct;

}

double Call_option::delta(double S, double t, double sigma, double r)
{
    double d1 = 1/(sigma * sqrt(T-t)) * ( log(S/K) + (r+ pow(sigma,2)/2)*(T-t) );
    return normalCDF(d1);

}

double Call_option::theta(double S, double t, double sigma, double r)
{
    double d1 = 1/(sigma * sqrt(T-t)) * ( log(S/K) + (r+ pow(sigma,2)/2)*(T-t) );
    double d2 = d1 - sigma*sqrt(T-t);
    double tht = -(1/sqrt(2*M_PI))*(S * exp(-pow(d1,2)/2) * sigma / (2*sqrt(T-t))) - r*K*exp(-r*(T-t))*normalCDF(d2);
    return tht;
}

double Call_option::gamma(double S, double t, double sigma, double r)
{
    double d1 = 1/(sigma * sqrt(T-t)) * ( log(S/K) + (r+ pow(sigma,2)/2)*(T-t) );
    double gam = exp(-pow(d1,2))/(S*sigma*sqrt(2*M_PI*(T-t)));
    return gam;
}

double Call_option::vega(double S, double t, double sigma, double r)
{
    double d1 = 1/(sigma * sqrt(T-t)) * ( log(S/K) + (r+ pow(sigma,2)/2)*(T-t) );
    double v = S*sqrt(T-t)*exp(-pow(d1,2)/2)/sqrt(2*M_PI);
    return v;
}

double Call_option::rho(double S, double t, double sigma, double r)
{
    double d1 = 1/(sigma * sqrt(T-t)) * ( log(S/K) + (r+ pow(sigma,2)/2)*(T-t) );
    double d2 = d1 - sigma*sqrt(T-t);
    double rh = (T-t)*K*exp(-r*(T-t))*normalCDF(d2);
    return rh;
}

double Call_option::impl_vol(double market_price, double sigma0, double S, double t, double r, double q, double eps)
{
    //use newton raphson
    double fun = BS(S, t, sigma0, r,q) - market_price;

    double deriv_fun = vega(S, t, sigma0, r);
    double h =  fun / deriv_fun;
    while (std::abs(h) >= eps)
    {
        fun = BS(S, t, sigma0, r,q) - market_price;
        deriv_fun = vega(S, t, sigma0, r);
        h = fun/deriv_fun;

        // x(i+1) = x(i) - f(x) / f'(x)
        sigma0 = sigma0 - h;
    }
    return sigma0;
}

Call_option::~Call_option()
{
    //dtor
}
