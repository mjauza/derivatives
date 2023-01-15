#include "American_option.h"
#include <string>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include "Tree_methods.h"

American_option::American_option(double T, double K, std::string type)
{
    if (type != "call" && type != "put"){
        throw std::invalid_argument( "type argument should be call/put" );
    }

    this -> T = T;
    this -> K = K;
    this -> type = type;
}

double American_option::call_tree_BS(int N, double r, double sigma, double S0)
{
    double delta_t = Tree_methods::get_delta_t(this -> T, N);
    double R = exp(r*delta_t);
    double u = Tree_methods::get_u(delta_t,sigma);
    double d = Tree_methods::get_d(delta_t, sigma);
    double p = Tree_methods::get_p(R,u,d);

    std::vector<std::vector<double>> lat = Tree_methods::get_lattice(delta_t, sigma, S0, N);

    double C_price[N+1][N+1];
    std::vector<double> ST_v = lat[lat.size()-1];
    for(int i=0; i<N+1; i++){
        double ST = ST_v[i];
        C_price[i][N] = std::max(ST - K, 0.0);
    }


    for(int i = N-1; i>=0; i--){
        for(int j=0; j<=i; j++){
            double c = p*C_price[j][i+1] + (1 - p)*C_price[j+1][i+1];
            double St = lat[i][j];
            double cur_c = std::max(St - K, 0.0);
            double exp_c = exp(-r*delta_t)*c;
            C_price[j][i] = std::max(exp_c, cur_c);
        }
    }

    return C_price[0][0];
}

double American_option::put_tree_BS(int N, double r, double sigma, double S0)
{
    double delta_t = Tree_methods::get_delta_t(this -> T, N);
    double R = exp(r*delta_t);
    double u = Tree_methods::get_u(delta_t,sigma);
    double d = Tree_methods::get_d(delta_t, sigma);
    double p = Tree_methods::get_p(R,u,d);

    std::vector<std::vector<double>> lat = Tree_methods::get_lattice(delta_t, sigma, S0, N);

    double C_price[N+1][N+1];
    std::vector<double> ST_v = lat[lat.size()-1];
    for(int i=0; i<N+1; i++){
        double ST = ST_v[i];
        C_price[i][N] = std::max(K - ST, 0.0);
    }


    for(int i = N-1; i>=0; i--){
        for(int j=0; j<=i; j++){
            double c = p*C_price[j][i+1] + (1 - p)*C_price[j+1][i+1];
            double St = lat[i][j];
            double cur_c = std::max(K - St, 0.0);
            double exp_c = exp(-r*delta_t)*c;
            C_price[j][i] = std::max(exp_c, cur_c);
        }
    }

    return C_price[0][0];
}

double American_option::tree_BS(int N, double r, double sigma, double S0)
{
    double price;
    if(this -> type == "call"){
        price = call_tree_BS(N,r,sigma,S0);
    }else{
        price = put_tree_BS(N,r,sigma,S0);
    }
    return price;
}

American_option::~American_option()
{
    //dtor
}
