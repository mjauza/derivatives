#include "Exchange_option.h"
#include <cmath>
Exchange_option::Exchange_option(double K, double T)
{
    this -> K = K;
    this -> T = T;
}

double Exchange_option::normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x/std::sqrt(2))/2;
}

double Exchange_option::sigma2_YX(double sigma_X, double sigma_Y, double rho)
{
    return pow(sigma_X,2) - 2*rho*sigma_X*sigma_Y + pow(sigma_Y,2);
}

double Exchange_option::dX(double X, double Y, double sigma_X, double sigma_Y, double rho, double t)
{
    double sigma2YX = sigma2_YX(sigma_X, sigma_Y, rho);
    double sigmaYX = sqrt(sigma2YX);
    return (log(X/Y) + sigma2YX*(T-t)) / (sigmaYX*sqrt(T-t));
}

double Exchange_option::dY(double X, double Y, double sigma_X, double sigma_Y, double rho, double t)
{
    double d_X = dX(X,Y,sigma_X,sigma_Y,rho,t);
    double sigmaYX = sqrt(sigma2_YX(sigma_X, sigma_Y, rho));
    return d_X - sigmaYX*sqrt(T-t);
}

double Exchange_option::BS(double X, double Y, double sigma_X, double sigma_Y, double rho, double t)
{
    double d_X = dX(X,Y,sigma_X,sigma_Y,rho,t);
    double d_Y = dY(X,Y,sigma_X,sigma_Y,rho,t);
    return X*normalCDF(d_X) - Y*normalCDF(d_Y);

}

Exchange_option::~Exchange_option()
{
    //dtor
}
