#ifndef CALL_OPTION_H
#define CALL_OPTION_H
#include <cmath>

class Call_option
{
    public:
        Call_option(double K, double T);
        virtual ~Call_option();
        double BS(double S, double t, double sigma, double r, double q = 0);
        double delta(double S, double t, double sigma, double r);
        double theta(double S, double t, double sigma, double r);
        double gamma(double S, double t, double sigma, double r);
        double vega(double S, double t, double sigma, double r);
        double rho(double S, double t, double sigma, double r);
        double impl_vol(double market_price, double sigma0, double S, double t, double r, double q = 0, double eps = pow(10, -8));
        double normalCDF(double x);

    protected:
        double K;
        double T;

    private:
};

#endif // CALL_OPTION_H
