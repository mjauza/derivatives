#include "Tree_methods.h"
#include <vector>
#include <cmath>

Tree_methods::Tree_methods()
{

}

std::vector<std::vector<double>> Tree_methods::get_lattice(double delta_t, double sigma, double S0, int N)
{
    double u = get_u(delta_t, sigma);
    double d = get_d(delta_t, sigma);

    std::vector<std::vector<double>> lattice;
    lattice.push_back({S0});
    for(int i=1; i<=N; i++){
        std::vector<double> price_t;
        std::vector<double> price_t1 = lattice[i-1];
        for(int j = 0; j<=i; j++){
            double price;
            if(j==i){
                price = price_t1[j-1]*d;
            }else{
                price = price_t1[j]*u;
            }
            price_t.push_back(price);
        }
        lattice.push_back(price_t);
    }

    return lattice;
}

double Tree_methods::get_delta_t(double T, int N)
{
    return T/N;
}

double Tree_methods::get_u(double delta_t, double sigma)
{
    return exp(sigma*sqrt(delta_t));
}

double Tree_methods::get_d(double delta_t, double sigma)
{
    double u = get_u(delta_t, sigma);
    return 1/u;
}

double Tree_methods::get_p(double R, double u, double d)
{
    return (R - d)/(u - d);
}


Tree_methods::~Tree_methods()
{
    //dtor
}
