#ifndef EXCHANGE_OPTION_H
#define EXCHANGE_OPTION_H


class Exchange_option
{
    public:
        //payout = max(X_T - Y_T, 0)
        Exchange_option(double K, double T);
        double normalCDF(double x);
        double BS(double X, double Y, double sigma_X, double sigma_Y, double rho, double t);
        virtual ~Exchange_option();

    protected:
        double K;
        double T;

    private:
        double sigma2_YX(double sigma_X, double sigma_Y, double rho);
        double dX(double X, double Y, double sigma_X, double sigma_Y, double rho, double t);
        double dY(double X, double Y, double sigma_X, double sigma_Y, double rho, double t);
};

#endif // EXCHANGE_OPTION_H
