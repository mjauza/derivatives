#ifndef LOOKBACK_OPTION_H
#define LOOKBACK_OPTION_H
#include <string>

class Lookback_option
{
    public:
        Lookback_option(double X, double T, std::string type);
        double normalCDF(double x);
        double BS(double S, double m_M, double sigma, double t, double r);
        virtual ~Lookback_option();

    protected:
        double X; //strike
        double T; //maturity
        std::string type;

    private:
        double H(double S, double K, double sigma, double tau, double r);
        double d(double S, double K, double sigma, double tau, double r);
        double BS_call(double S, double M, double sigma, double t, double r);
        double h(double S, double K, double sigma, double tau, double r);
        double BS_put(double S, double m, double t, double sigma, double r);
};

#endif // LOOKBACK_OPTION_H
