#include "epcsaft_core_internal.h"

// EqID: h_res
double hres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    double Z = Z_cpp(t, rho, x, cppargs);
    double dares_dt = dadt_cpp(t, rho, std::move(x), cppargs);
    return (-t * dares_dt + (Z - 1.0)) * kb * N_AV * t;
}
