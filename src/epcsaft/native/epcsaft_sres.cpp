#include "epcsaft_core_internal.h"

// EqID: s_res
double sres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    double gres = gres_cpp(t, rho, x, cppargs);
    double hres = hres_cpp(t, rho, std::move(x), cppargs);
    return (hres - gres) / t;
}
