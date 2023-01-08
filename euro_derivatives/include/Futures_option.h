#ifndef FUTURES_OPTION_H
#define FUTURES_OPTION_H
#include <string>

class Futures_option
{
    public:
        Futures_option(double K, double T, std::string type);
        double normalCDF(double x);
        virtual ~Futures_option();
        double BS(double f, double sigma, double t, double r);

    protected:
        double K;
        double T;
        std::string type;
        double put_BS(double f, double sigma, double t, double r);
        double call_BS(double f, double sigma, double t, double r);

    private:
        double d1(double f, double sigma, double t);
        double d2(double f, double sigma, double t);
};

#endif // FUTURES_OPTION_H
