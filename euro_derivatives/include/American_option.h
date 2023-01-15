#ifndef AMERICAN_OPTION_H
#define AMERICAN_OPTION_H
#include <string>

class American_option
{
    public:
        American_option(double T, double K, std::string type);
        double tree_BS(int N, double r, double sigma, double S0);
        virtual ~American_option();

    protected:
        double T;
        double K;
        std::string type;
        double call_tree_BS(int N, double r, double sigma, double S0);
        double put_tree_BS(int N, double r, double sigma, double S0);

    private:
};

#endif // AMERICAN_OPTION_H
