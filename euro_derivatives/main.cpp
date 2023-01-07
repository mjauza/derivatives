#include <iostream>
#include "Call_option.h"

using namespace std;

int main()
{
    cout << "Hello world!" << endl;
    double K = 90;
    double T = 2;
    Call_option opt = Call_option(K, T);

    double S = 100;
    double t = 0;
    double sigma = 0.1;
    double r = 0.01;
    double bs_price = opt.BS(S,t,sigma,r);

    std::cout << "BS price: " << bs_price << "\n";
    return 0;
}
