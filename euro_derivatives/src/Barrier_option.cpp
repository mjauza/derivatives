#include "Barrier_option.h"
#include <string>
#include <cmath>
#include "Call_option.h"
#include <stdexcept>
#include <iostream>
Barrier_option::Barrier_option(double K, double T, std::string type, double B, std::string barrier_type)
{

    if (type != "call"){
        throw std::invalid_argument( "type argument should be call. Since at the moment only call options are implemnted" );
    }

    if( !((type == "call") || (type == "put")) ){
        throw std::invalid_argument( "type argument should be either call/put" );
    }

    if( !((barrier_type == "do") || (barrier_type == "di") || (barrier_type == "uo") || (barrier_type == "ui")) ){
        throw std::invalid_argument( "barrier_type argument should be either do/di/uo/ui" );
    }
    this -> K = K;
    this -> T = T;
    this -> B = B;
    this -> type = type;
    this -> barrier_type = barrier_type;
}

double Barrier_option::normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x/std::sqrt(2))/2;
}

double Barrier_option::d1(double S, double r, double sigma, double t)
{
    double K_num = std::max(B,K);
    return (log(S/K_num) + (r + pow(sigma,2)/2)*(T-t)) / (sigma * sqrt(T-t));
}

double Barrier_option::d2(double S, double r, double sigma, double t)
{
    return d1(S,r,sigma,t) - sigma*sqrt(T-t);
}

double Barrier_option::d3(double S, double r, double sigma, double t)
{
    return d1(S,r,sigma,t) + 2/(sigma*sqrt(T-t))*log(B/S);
}

double Barrier_option::d4(double S, double r, double sigma, double t)
{
    return d2(S,r,sigma,t) + 2/(sigma*sqrt(T-t))*log(B/S);
}

double Barrier_option::BS_call_do(double S, double r, double sigma, double t)
{
    double delta_num = 2*r/pow(sigma,2);
    double d_1 = d1(S,r,sigma,t);
    double d_3 = d3(S,r,sigma,t);
    double d_2 = d2(S,r,sigma,t);
    double d_4 = d4(S,r,sigma,t);
    return S*(normalCDF(d_1) - pow(B/S, delta_num+1)*normalCDF(d_3)) - K*exp(-r*(T-t))*(normalCDF(d_2) - pow(B/S,delta_num-1)*normalCDF(d_4));
}

double Barrier_option::BS_call_di(double S, double r, double sigma, double t)
{
    double price = pow(B,2)/S;
    double delta_num = 2*r/pow(sigma,2);
    Call_option call_euro = Call_option(this->K, T-t);
    double call_price = call_euro.BS(price, 0,sigma,r);
    return pow(B/S, delta_num-1)*call_price;

}

double Barrier_option::BS_call_uo(double S, double r, double sigma, double t)
{
    double delta_num = 2*r/pow(sigma,2);

    Call_option euro_call_a = Call_option(K,T-t);
    double price_a_1 = euro_call_a.BS(S, 0, sigma, r);
    double price_a_2 = euro_call_a.BS(pow(B,2)/S, 0, sigma, r);

    Call_option euro_call_b = Call_option(B,T-t);
    double price_b_1 = euro_call_b.BS(S, 0, sigma, r);
    double price_b_2 = euro_call_b.BS(pow(B,2)/S, 0, sigma, r);

    return (price_a_1 - pow(B/S, delta_num - 1)*price_a_2) - (price_b_1 - pow(B/S, delta_num - 1)*price_b_2);
}

double Barrier_option::BS_call_ui(double S, double r, double sigma, double t)
{
    double delta_num = 2*r/pow(sigma,2);
    Call_option euro_call_1 = Call_option(K,T-t);
    double euro_call_1_price = euro_call_1.BS(pow(B,2)/S, 0, sigma, r);
    Call_option euro_call_2 = Call_option(B,T-t);
    double euro_call_2_price = euro_call_2.BS(S, 0, sigma, r);
    Call_option euro_call_3 = Call_option(B,T-t);
    double euro_call_3_price = euro_call_3.BS(pow(B,2)/S, 0, sigma, r);
    return pow(B/S, delta_num-1)*euro_call_1_price + euro_call_2_price - pow(B/S,delta_num-1)*euro_call_3_price;
}

double Barrier_option::BS(double S, double r, double sigma, double t)
{
    double price;
    if(barrier_type == "do"){
        price = BS_call_do(S, r, sigma, t);
    }else if(barrier_type == "di"){
        price = BS_call_di(S, r, sigma, t);
    } else if (barrier_type == "uo"){
        price = BS_call_uo(S, r, sigma, t);
    } else {
        price = BS_call_ui(S, r, sigma, t);
    }
    return price;
}

Barrier_option::~Barrier_option()
{
    //dtor
}
