#ifndef BARRIER_OPTION_H
#define BARRIER_OPTION_H
#include <string>
#include "Call_option.h"
class Barrier_option
{
    public:
        Barrier_option(double K, double T, std::string type, double B, std::string barrier_type);
        double normalCDF(double x);
        double BS(double S, double r, double sigma, double t);
        virtual ~Barrier_option();

    protected:
        double K;
        double T;
        std::string type;
        double B;
        std::string barrier_type;

    private:
        double BS_call_do(double S, double r, double sigma, double t);
        double BS_call_di(double S, double r, double sigma, double t);
        double BS_call_uo(double S, double r, double sigma, double t);
        double BS_call_ui(double S, double r, double sigma, double t);

        double d1(double S, double r, double sigma, double t);
        double d2(double S, double r, double sigma, double t);
        double d3(double S, double r, double sigma, double t);
        double d4(double S, double r, double sigma, double t);

};

#endif // BARRIER_OPTION_H
