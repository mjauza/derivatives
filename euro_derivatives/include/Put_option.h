#ifndef PUT_OPTION_H
#define PUT_OPTION_H


class Put_option
{
    public:
        Put_option(double K, double T);
        double BS(double S, double t, double sigma, double r, double q = 0);
        double delta(double S, double t, double sigma, double r);
        double theta(double S, double t, double sigma, double r);
        double gamma(double S, double t, double sigma, double r);
        double vega(double S, double t, double sigma, double r);
        double rho(double S, double t, double sigma, double r);
        double normalCDF(double x);
        virtual ~Put_option();

    protected:
        double K;
        double T;

    private:
};

#endif // PUT_OPTION_H
