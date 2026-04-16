#include "epcsaft_contrib_internal.h"

using namespace thermo_detail;

DispersionPolynomialState dispersion_polynomials_cpp(double m_avg, double eta) {
    DispersionPolynomialState state;
    double c1 = (m_avg - 1.0) / m_avg;
    double c2 = (m_avg - 2.0) / m_avg;
    for (size_t i = 0; i < state.a.size(); ++i) {
        state.a[i] = kDispersionA0[i] + c1 * kDispersionA1[i] + c1 * c2 * kDispersionA2[i];
        state.b[i] = kDispersionB0[i] + c1 * kDispersionB1[i] + c1 * c2 * kDispersionB2[i];
        state.I1 += state.a[i] * std::pow(eta, static_cast<int>(i));
        state.I2 += state.b[i] * std::pow(eta, static_cast<int>(i));
        if (i > 0) {
            state.dI1_deta += state.a[i] * static_cast<double>(i) * std::pow(eta, static_cast<int>(i - 1));
            state.dI2_deta += state.b[i] * static_cast<double>(i) * std::pow(eta, static_cast<int>(i - 1));
        }
        state.dEtaI1_deta += state.a[i] * static_cast<double>(i + 1) * std::pow(eta, static_cast<int>(i));
        state.dEtaI2_deta += state.b[i] * static_cast<double>(i + 1) * std::pow(eta, static_cast<int>(i));
    }
    state.C1 = 1.0 / (1.0
        + m_avg * (8.0 * eta - 2.0 * eta * eta) / std::pow(1.0 - eta, 4.0)
        + (1.0 - m_avg) * (20.0 * eta - 27.0 * eta * eta + 12.0 * std::pow(eta, 3.0) - 2.0 * std::pow(eta, 4.0))
            / std::pow((1.0 - eta) * (2.0 - eta), 2.0));
    state.C2 = -state.C1 * state.C1 * (
        m_avg * (-4.0 * eta * eta + 20.0 * eta + 8.0) / std::pow(1.0 - eta, 5.0)
        + (1.0 - m_avg) * (2.0 * std::pow(eta, 3.0) + 12.0 * eta * eta - 48.0 * eta + 40.0)
            / std::pow((1.0 - eta) * (2.0 - eta), 3.0));
    return state;
}
