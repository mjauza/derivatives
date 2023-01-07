#ifndef CALL_OPTION_H
#define CALL_OPTION_H


class Call_option
{
    public:
        Call_option(double K, double T);
        virtual ~Call_option();
        double BS(double S, double t, double sigma, double r);
        double normalCDF(double x);

    protected:
        double K;
        double T;

    private:
};

#endif // CALL_OPTION_H
