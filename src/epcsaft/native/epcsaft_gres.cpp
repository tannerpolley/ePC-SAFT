#include "epcsaft_core_internal.h"

// EqID: g_res_from_ares
double gres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    double ares = ares_cpp(t, rho, x, cppargs);
    double Z = Z_cpp(t, rho, std::move(x), cppargs);
    return (ares + (Z - 1.0) - std::log(Z)) * kb * N_AV * t;
}
