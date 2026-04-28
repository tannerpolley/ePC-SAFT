#include "epcsaft_core_internal.h"

// EqID: h_res
double hres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    double Z = Z_cpp(t, rho, x, cppargs);
    double dares_dt = dadt_cpp(t, rho, std::move(x), cppargs);
    return (-t * dares_dt + (Z - 1.0)) * kb * N_AV * t;
}

// EqID: g_res_from_ares
double gres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    double ares = ares_cpp(t, rho, x, cppargs);
    double Z = Z_cpp(t, rho, std::move(x), cppargs);
    return (ares + (Z - 1.0) - std::log(Z)) * kb * N_AV * t;
}

// EqID: s_res
double sres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    double gres = gres_cpp(t, rho, x, cppargs);
    double hres = hres_cpp(t, rho, std::move(x), cppargs);
    return (hres - gres) / t;
}
