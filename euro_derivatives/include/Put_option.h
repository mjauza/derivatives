#ifndef PUT_OPTION_H
#define PUT_OPTION_H
#include <cmath>
class Put_option
{
    public:
        Put_option(double K, double T);
        double BS(double S, double t, double sigma, double r, double q = 0);
        double MC_BS(double S, double t, double sigma, double r, int N);
        double tree_BS(double sigma, double S0, double r, int N);

        double delta(double S, double t, double sigma, double r);
        double theta(double S, double t, double sigma, double r);
        double gamma(double S, double t, double sigma, double r);
        double vega(double S, double t, double sigma, double r);
        double rho(double S, double t, double sigma, double r);
        double impl_vol(double market_price, double sigma0, double S, double t, double r, double q = 0, double eps = pow(10, -8));
        double normalCDF(double x);
        virtual ~Put_option();

    protected:
        double K;
        double T;

    private:
};

#endif // PUT_OPTION_H
