#include "Lookback_option.h"
#include <string>
#include <iostream>
#include <stdexcept>
#include <cmath>

Lookback_option::Lookback_option(double X, double T, std::string type)
{
    if ( !((type == "call") || (type == "put")) ){
        throw std::invalid_argument( "type argument should be call/put" );
    }

    this -> X = X;
    this -> T = T;
    this -> type = type;
}

double Lookback_option::normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x/std::sqrt(2))/2;
}

double Lookback_option::d(double S, double K, double sigma, double tau, double r)
{
    return (log(S/K) + (r + pow(sigma,2)/2)*tau) / (sigma * sqrt(tau));
}
double Lookback_option::H(double S, double K, double sigma, double tau, double r)
{
    double d_num = d(S,K,sigma,tau,r);
    double a_num = S*normalCDF(d_num) - exp(-r*tau)*K*normalCDF(d_num - sigma*sqrt(tau));
    double b_num = exp(-r*tau)*pow(sigma,2)/(2*r)*S*(exp(r*tau)*normalCDF(d_num) - pow(S/K, -2*r/pow(sigma,2))*normalCDF(d_num - 2*r/sigma*sqrt(tau)));
    return a_num + b_num;
}

double Lookback_option::BS_call(double S, double M, double sigma, double t, double r)
{
    double a_num = exp(-r*(T-t))*std::max(M - X, 0.0);
    double b_num = H(S, std::max(M,X), sigma, T-t, r);
    return a_num+b_num;
}

double Lookback_option::h(double S, double K, double sigma, double tau, double r)
{
    double d_num = d(S,K,sigma,tau,r);
    double a_num = exp(-r*tau)*K*normalCDF(-d_num + sigma*sqrt(tau)) - S*normalCDF(-d_num);
    double b_num = exp(-r*tau)*pow(sigma, 2)/(2*r)*S;
    double c_num = pow(S/K, -2*r/pow(sigma,2))*normalCDF(-d_num + 2*r/sigma*sqrt(tau)) - exp(r*tau)*normalCDF(-d_num);
    return a_num + b_num*c_num;
}

double Lookback_option::BS_put(double S, double m, double t, double sigma, double r)
{
    double a_num = exp(-r*(T-t))*std::max(X-m,0.0);
    double b_num = h(S,std::min(m,X), sigma, T-t,r);
    return a_num + b_num;
}

double Lookback_option::BS(double S, double m_M, double sigma, double t, double r)
{
    double price;
    if(this -> type == "call"){
        price = BS_call(S,m_M,sigma,t,r);
    }else{
        price = BS_put(S,m_M,t,sigma,r);
    }
    return price;

}

Lookback_option::~Lookback_option()
{
    //dtor
}
