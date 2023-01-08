#ifndef CHOOSER_OPTION_H
#define CHOOSER_OPTION_H
#include "Call_option.h"

class Chooser_option
{
    public:
        Chooser_option(double K, double T);
        double BS_Tc(double S_Tc, double Tc, double sigma, double r, double q = 0);
        virtual ~Chooser_option();

    protected:
        double K;
        double T;

    private:
        Call_option c;
};

#endif // CHOOSER_OPTION_H
