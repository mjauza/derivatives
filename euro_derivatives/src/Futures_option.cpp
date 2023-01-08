#include <iostream>
#include "Futures_option.h"
#include <stdexcept>
#include <cmath>

Futures_option::Futures_option(double K, double T, std::string type)
{
    this -> K = K;
    this -> T = T;

    if( !((type == "call") || (type == "put")) ){
        throw std::invalid_argument( "type argument should be either call or put" );
    }
    this -> type = type;
}

double Futures_option::normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x/std::sqrt(2))/2;
}

double Futures_option::d1(double f, double sigma, double t)
{
    return (log(f/K) + pow(sigma,2)*(T-t)) / (sigma * sqrt(T-t));
}

double Futures_option::d2(double f, double sigma, double t)
{
    double d_1 = this -> d1(f,sigma,t);
    return d_1 - sigma*sqrt(T-t);
}

double Futures_option::put_BS(double f, double sigma, double t, double r)
{
    double d_2 = d2(f,sigma,t);
    double d_1 = d1(f,sigma,t);
    double p = exp(-r*(T-t))*(K*normalCDF(-d_2) - f*normalCDF(-d_1));
    return p;
}

double Futures_option::call_BS(double f, double sigma, double t, double r)
{
    double d_2 = d2(f,sigma,t);
    double d_1 = d1(f,sigma,t);
    double c = exp(-r*(T-t))*(f*normalCDF(d_1) - K*normalCDF(d_2));
    return c;
}

double Futures_option::BS(double f, double sigma, double t, double r)
{
    if(this->type == "call"){
        return call_BS(f, sigma, t, r);
    } else {
        return put_BS(f, sigma, t, r);
    }
}

Futures_option::~Futures_option()
{
    //dtor
}
