#include "Simulation_methods.h"
#include <random>
#include <cmath>
#include <string>
#include <map>
#include <tuple>
#include <vector>
Simulation_methods::Simulation_methods()
{
    //ctor
}

std::tuple<std::vector<double>, std::vector<double>> Simulation_methods::BM(int N, double T)
{
    double t_delta = T/N;
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,sqrt(t_delta));
    std::vector<double> path_values;
    std::vector<double> path_times;
    path_values.push_back(0);
    path_times.push_back(0);
    for(int i=1; i<=N; i++){
        double inc = distribution(generator);
        path_values.push_back(path_values[i-1] + inc);
        path_times.push_back(path_times[i-1] + t_delta);
    }
    return {path_times, path_values};
}

std::tuple<std::vector<double>, std::vector<double>> Simulation_methods::geom_BM(int N, double T, double mu, double sigma, double S0)
{
    auto [bm_times, bm_values] = BM(N,T);
    std::vector<double> path_values;
    double drift = mu - pow(sigma,2)/2;
    for(int i= 0; i<bm_times.size(); i++){
        double Wt = bm_values[i];
        double t = bm_times[i];
        double St = S0*exp(drift*t + sigma*Wt);
        path_values.push_back(St);
    }
    return {bm_times, path_values};

}

Simulation_methods::~Simulation_methods()
{
    //dtor
}
