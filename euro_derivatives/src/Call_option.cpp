#include "Call_option.h"
#include <cmath>
#define _USE_MATH_DEFINES
#include <iostream>
#include "Simulation_methods.h"
#include <numeric>
#include "Tree_methods.h"
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

    // avergae
    double C_sum = std::accumulate(CT.begin(), CT.end(), 0);
    double C_avg = C_sum / N;

    // discount
    double Ct = exp(-r*(T-t))*C_avg;
    return Ct;

}

double Call_option::tree_BS(double sigma, double S0, double r, int N)
{
    double delta_t = Tree_methods::get_delta_t(this -> T, N);
    double R = exp(r*delta_t);
    double u = Tree_methods::get_u(delta_t,sigma);
    double d = Tree_methods::get_d(delta_t, sigma);
    double p = Tree_methods::get_p(R,u,d);

    std::vector<std::vector<double>> lat = Tree_methods::get_lattice(delta_t, sigma, S0, N);

    // get last layer and calculate price
    std::vector<double> CT;
    for(double s : lat[lat.size()-1]){
        double c = std::max(s - this->K, 0.0);
        CT.push_back(c);
    }

    double C_price[N+1][N+1];
    for(int i=0; i<N+1; i++){
        C_price[i][N] = CT[i];
    }

    // create a latice of call prices
    for(int i = N-1; i>=0; i--){
        for(int j=0; j<=i; j++){
            double c = p*C_price[j][i+1] + (1 - p)*C_price[j+1][i+1];
            C_price[j][i] = exp(-r*delta_t)*c;
        }
    }

    double final_price = C_price[0][0];
    //std::cout << final_price;

    return final_price;

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
