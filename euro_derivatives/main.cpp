#include <iostream>
#include "Call_option.h"
#include "Put_option.h"
#include "Futures_option.h"
#include "Chooser_option.h"
#include "Exchange_option.h"
#include "Barrier_option.h"
#include "Lookback_option.h"
#include "Simulation_methods.h"
#include <tuple>
#include <vector>
#include "Tree_methods.h"

using namespace std;
void test_BM();
void test_geom_BM();
void test_EU_call();
void test_EU_put();
void test_lattice();

int main()
{
    //test_BM();
    //test_geom_BM();

    test_EU_call();
    //test_EU_put();

    //test_lattice();

    /*
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

    // barrier options
    Barrier_option bar_call_do = Barrier_option(K, T, "call", 70, "do");
    double bar_call_do_price = bar_call_do.BS(80, r,sigma,t);
    Barrier_option bar_call_di = Barrier_option(K, T, "call", 70, "di");
    double bar_call_di_price = bar_call_di.BS(80, r,sigma,t);
    Barrier_option bar_call_uo = Barrier_option(90, T, "call", 100, "uo");
    double bar_call_uo_price = bar_call_uo.BS(80, r,sigma,t);
    Barrier_option bar_call_ui = Barrier_option(K, T, "call", 95, "ui");
    double bar_call_ui_price = bar_call_ui.BS(80, r,sigma,t);

    std::cout << "bar_call_do_price = " << bar_call_do_price << "\n";
    std::cout << "bar_call_di_price = " << bar_call_di_price << "\n";
    std::cout << "bar_call_uo_price = " << bar_call_uo_price << "\n";
    std::cout << "bar_call_ui_price = " << bar_call_ui_price << "\n";

    // lookback option
    double S_lb = 100;
    double m_M = S_lb;
    double X_lb = 90;
    Lookback_option lb_call = Lookback_option(X_lb,3,"call");
    double lb_call_price = lb_call.BS(S_lb,m_M,sigma,0,r);
    Lookback_option lb_put = Lookback_option(X_lb,3,"put");
    double lb_put_price = lb_put.BS(S_lb,m_M,sigma,0,r);

    std::cout <<"Lookback call price: " << lb_call_price << "\n";
    std::cout <<"Lookback put price: " << lb_put_price << "\n";

    */

    ///////////////////////////////////////
    return 0;
}

void test_BM()
{
    //std::tuple<std::vector<double>, std::vector<double>> bm_path_tuple = Simulation_methods::BM(100, 1.0);
    //std::vector<double> bm_times = std::get<0>(bm_path_tuple);
    //std::vector<double> bm_values = std::get<1>(bm_path_tuple);
    auto [bm_times, bm_values] = Simulation_methods::BM(100, 1.0);
    for(int i = 0; i < bm_times.size(); i++){
        std::cout << "BM time = " << bm_times[i] << " ; value = " << bm_values[i] << "\n";
    }
}
void test_geom_BM()
{
    auto [bm_times, bm_values] = Simulation_methods::geom_BM(100, 1.0,0.1,0.2,100);
    for(int i = 0; i < bm_times.size(); i++){
        std::cout << "Geometric BM time = " << bm_times[i] << " ; value = " << bm_values[i] << "\n";
    }
}

void test_EU_call()
{
    double S = 100;
    double K = 70;
    double T = 2;
    double r = 0.02;
    double sigma = 0.02;
    const double t = 0;
    double N = 10000;
    Call_option opt = Call_option(K, T);
    double bs_price = opt.BS(S,t,sigma,r);
    double mc_price = opt.MC_BS(S,t,sigma,r,N);
    double tree_price = opt.tree_BS(sigma, S, r, 10);
    std::cout << "BS price = " << bs_price << "\n";
    std::cout << "MC price = " << mc_price << "\n";
    std::cout << "Tree price = " << tree_price << "\n";

}

void test_EU_put()
{
    double S = 65;
    double K = 70;
    double T = 2;
    double r = 0.02;
    double sigma = 0.02;
    double t = 0;
    int N = 100000;
    Put_option opt = Put_option(K, T);
    double bs_price = opt.BS(S,t,sigma,r);
    double mc_price = opt.MC_BS(S,t,sigma,r,N);
    double tree_price = opt.tree_BS(sigma, S, r, 100);
    std::cout << "BS price = " << bs_price << "\n";
    std::cout << "MC price = " << mc_price << "\n";
    std::cout << "Tree price = " << tree_price << "\n";
}

void test_lattice()
{
    std::vector<std::vector<double>> lat = Tree_methods::get_lattice( 0.1, 0.2, 100, 10);
    for(int t = 0; t<lat.size(); t++){
        std::vector<double> price_v = lat[t];
        for(int i = 0; i<price_v.size(); i++){
            std::cout << "t = " << t << " i = " << i  << " price = " << price_v[i] << "\n";
        }
    }
}
