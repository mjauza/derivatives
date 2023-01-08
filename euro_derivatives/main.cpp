#include <iostream>
#include "Call_option.h"
#include "Put_option.h"
#include "Futures_option.h"
#include "Chooser_option.h"
#include "Exchange_option.h"

using namespace std;

int main()
{


    double K = 90;
    double T = 2;
    Call_option call_opt = Call_option(K, T);
    Put_option put_opt = Put_option(K, T);

    double S = 100;
    double t = 0;
    double Tc = 1;
    double S_Tc = 95;
    double sigma = 0.1;
    double r = 0.01;
    double call_bs_price = call_opt.BS(S,t,sigma,r);
    double put_bs_price = put_opt.BS(S,t,sigma,r);

    std::cout << "Call BS price: " << call_bs_price << "\n";
    std::cout << "Put BS price: " << put_bs_price << "\n";

    // get implied vol
    double call_imp_vol = call_opt.impl_vol(13, 0.4, S, t, r);
    double put_imp_vol = put_opt.impl_vol(1.3, 0.4, S, t, r);
    std::cout << "Call BS implied vol: " << call_imp_vol << "\n";
    std::cout << "Put BS implied vol: " << put_imp_vol << "\n";

    // get deltas
    double call_delta = call_opt.delta(S, t, sigma, r);
    double put_delta = put_opt.delta(S, t, sigma, r);

    std::cout << "Call delta: " << call_delta << "\n";
    std::cout << "Put delta: " << put_delta << "\n";

    // get theta
    double call_theta = call_opt.theta(S, t, sigma, r);
    double put_theta = put_opt.theta(S, t, sigma, r);
    std::cout << "Call theta: " << call_theta << "\n";
    std::cout << "Put theta: " << put_theta << "\n";

    // get gamma
    double call_gamma = call_opt.gamma(S, t, sigma, r);
    double put_gamma = put_opt.gamma(S, t, sigma, r);
    std::cout << "Call gamma: " << call_gamma << "\n";
    std::cout << "Put gamma: " << put_gamma << "\n";

    // get vega
    double call_vega = call_opt.vega(S, t, sigma, r);
    double put_vega = put_opt.vega(S, t, sigma, r);
    std::cout << "Call vega: " << call_vega << "\n";
    std::cout << "Put vega: " << put_vega << "\n";

    // get rho
    double call_rho = call_opt.rho(S, t, sigma, r);
    double put_rho = put_opt.rho(S, t, sigma, r);
    std::cout << "Call rho: " << call_rho << "\n";
    std::cout << "Put rho: " << put_rho << "\n";

    // futures options
    Futures_option call_future_opt = Futures_option(K, T, "call");
    double f = 100;
    double call_future_bs = call_future_opt.BS(f,sigma,t,r);
    std::cout << "Future call BS: " << call_future_bs << "\n";

    // chooser option
    Chooser_option choos_opt = Chooser_option(K,T);
    double choos_opt_bs = choos_opt.BS_Tc(S_Tc, Tc, sigma, r);
    std::cout << "Chooser option BS at Tc: " << choos_opt_bs << "\n";

    //exchange options
    Exchange_option ex_option = Exchange_option(K,T);
    double X = 100;
    double Y = 90;
    double sigma_X = 0.2;
    double sigma_Y = 0.2;
    double rho = 0.4;
    double ex_option_price = ex_option.BS(X, Y, sigma_X, sigma_Y, rho, t);
    std::cout << "Exchange option BS: " << ex_option_price << "\n";


    ///////////////////////////////////////
    return 0;
}
