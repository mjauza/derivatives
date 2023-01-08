#include <iostream>
#include "Call_option.h"
#include "Put_option.h"


using namespace std;

int main()
{
    cout << "Hello world!" << endl;
    double K = 90;
    double T = 2;
    Call_option call_opt = Call_option(K, T);
    Put_option put_opt = Put_option(K, T);

    double S = 100;
    double t = 0;
    double sigma = 0.1;
    double r = 0.01;
    double call_bs_price = call_opt.BS(S,t,sigma,r);
    double put_bs_price = put_opt.BS(S,t,sigma,r);

    std::cout << "Call BS price: " << call_bs_price << "\n";
    std::cout << "Put BS price: " << put_bs_price << "\n";

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
    return 0;
}
