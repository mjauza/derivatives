#ifndef TREE_METHODS_H
#define TREE_METHODS_H
#include <vector>

class Tree_methods
{
    public:
        Tree_methods();
        virtual ~Tree_methods();
        static std::vector<std::vector<double>> get_lattice(double delta_t, double sigma, double S0, int N);
        static double get_delta_t(double T, int N);
        static double get_u(double delta_t, double sigma);
        static double get_d(double delta_t, double sigma);
        static double get_p(double R,double u,double d);



    protected:





    private:
};

#endif // TREE_METHODS_H
