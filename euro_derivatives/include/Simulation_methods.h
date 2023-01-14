#ifndef SIMULATION_METHODS_H
#define SIMULATION_METHODS_H
#include <tuple>
#include <vector>

class Simulation_methods
{
    public:
        Simulation_methods();

        // standard BM
        static std::tuple<std::vector<double>, std::vector<double>> BM(int N, double T);
        static std::vector<double> BM_t(int N, double t);

        // geometric BM
        static std::tuple<std::vector<double>, std::vector<double>> geom_BM(int N, double T, double mu, double sigma, double S0);
        static std::vector<double> geom_BM_t(int N, double t, double mu, double sigma, double S0);



        // log normal dist
        virtual ~Simulation_methods();

    protected:

    private:
};

#endif // SIMULATION_METHODS_H
