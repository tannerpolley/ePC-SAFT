#include "epcsaft_core_internal.h"

using thermo_detail::DensityBracket;
using thermo_detail::DensityRootCandidate;
using thermo_detail::DensityScanPoint;
using thermo_detail::parameter_setup_detail::ion_diameter_cpp;

vector<double> density_scan_grid_cpp() {
    const double nu_min = 1e-13;
    const double nu_log_max = 5e-3;
    const double nu_linear_start = 5e-3;
    const double nu_max = 0.7405 - 1e-4;
    const int n_log = 24;
    const int n_linear = 256;

    vector<double> grid;
    grid.reserve(1 + n_log + n_linear);
    grid.push_back(nu_min);

    if (n_log > 0) {
        const double log_min = std::log(nu_min);
        const double log_max = std::log(nu_log_max);
        for (int i = 0; i < n_log; i++) {
            double frac = static_cast<double>(i) / static_cast<double>(n_log - 1);
            double nu = std::exp(log_min + frac * (log_max - log_min));
            if (nu > grid.back()) {
                grid.push_back(nu);
            }
        }
    }

    if (n_linear > 0) {
        for (int i = 0; i < n_linear; i++) {
            double frac = static_cast<double>(i) / static_cast<double>(n_linear - 1);
            double nu = nu_linear_start + frac * (nu_max - nu_linear_start);
            if (nu > grid.back()) {
                grid.push_back(nu);
            }
        }
    }

    return grid;
}

DensityScanPoint density_scan_point_cpp(double nu, double t, int ncomp, const vector<double> &x, double p, const add_args &cppargs) {
    DensityScanPoint point;
    point.nu = nu;
    point.rho = ::reduced_density_to_molar(nu, t, ncomp, x, cppargs);
    try {
        point.resid = ::density_root_residual_cpp(point.rho, t, p, x, cppargs);
        point.finite = std::isfinite(point.resid);
    }
    catch (const std::exception&) {
        point.resid = 0.0;
        point.finite = false;
    }
    return point;
}

std::string density_failure_message_cpp(const std::string &outcome, double t, double p, const vector<double> &x, int phase) {
    std::string phase_name = (phase == 1) ? "vapor" : "liquid";
    return outcome
        + " (phase=" + phase_name
        + ", T=" + std::to_string(t)
        + ", P=" + std::to_string(p)
        + ", ncomp=" + std::to_string(x.size())
        + ").";
}

vector<DensityBracket> density_brackets_cpp(const vector<DensityScanPoint> &points) {
    vector<DensityBracket> brackets;
    if (points.size() < 2) {
        return brackets;
    }

    for (size_t i = 1; i < points.size(); i++) {
        const DensityScanPoint &lo = points[i - 1];
        const DensityScanPoint &hi = points[i];
        if (!lo.finite || !hi.finite) {
            continue;
        }
        if (lo.resid * hi.resid < 0.0) {
            DensityBracket bracket;
            bracket.nu_lo = lo.nu;
            bracket.nu_hi = hi.nu;
            bracket.resid_lo = lo.resid;
            bracket.resid_hi = hi.resid;
            brackets.push_back(bracket);
        }
    }

    return brackets;
}

void refine_density_brackets_cpp(
    const DensityBracket &coarse,
    double t,
    int ncomp,
    const vector<double> &x,
    double p,
    const add_args &cppargs,
    vector<DensityBracket> &refined_brackets
) {
    const int refine_segments = 256;
    const double discontinuity_threshold = 1e12;

    vector<DensityScanPoint> refined_points;
    refined_points.reserve(refine_segments + 1);
    for (int i = 0; i <= refine_segments; i++) {
        double frac = static_cast<double>(i) / static_cast<double>(refine_segments);
        double nu = coarse.nu_lo + frac * (coarse.nu_hi - coarse.nu_lo);
        refined_points.push_back(density_scan_point_cpp(nu, t, ncomp, x, p, cppargs));
    }

    vector<DensityBracket> local_brackets = density_brackets_cpp(refined_points);
    for (const DensityBracket &candidate : local_brackets) {
        double nu_mid = 0.5 * (candidate.nu_lo + candidate.nu_hi);
        DensityScanPoint mid = density_scan_point_cpp(nu_mid, t, ncomp, x, p, cppargs);
        if (!mid.finite) {
            continue;
        }

        double min_abs = std::min(std::abs(candidate.resid_lo), std::abs(candidate.resid_hi));
        min_abs = std::min(min_abs, std::abs(mid.resid));
        if (min_abs > discontinuity_threshold) {
            continue;
        }

        refined_brackets.push_back(candidate);
    }
}

bool density_root_valid_cpp(
    double t,
    double p,
    const vector<double> &x,
    const add_args &cppargs,
    double rho,
    DensityRootCandidate *candidate
) {
    const double rel_tol = 1e-6;
    const double abs_tol = 1e-7;
    const double ultra_low_pressure_rel_tol = 2.5e-4;
    const double ultra_low_pressure_cutoff = 1e-3;

    double p_calc = 0.0;
    try {
        p_calc = p_cpp(t, rho, x, cppargs);
    }
    catch (const std::exception&) {
        return false;
    }
    if (!std::isfinite(p_calc)) {
        return false;
    }

    double abs_p_error = std::abs(p_calc - p);
    double rel_resid = abs_p_error / std::max(std::abs(p), 1e-300);
    bool pressure_ok = (rel_resid <= rel_tol || abs_p_error <= abs_tol);
    if (!pressure_ok && p <= ultra_low_pressure_cutoff) {
        pressure_ok = (rel_resid <= ultra_low_pressure_rel_tol);
    }
    if (!pressure_ok) {
        return false;
    }

    double h = std::max(1e-12, std::abs(rho) * 1e-6);
    double dpdrho = 0.0;
    try {
        double p_plus = p_cpp(t, rho + h, x, cppargs);
        if (!std::isfinite(p_plus)) {
            return false;
        }
        if (rho > h) {
            double p_minus = p_cpp(t, rho - h, x, cppargs);
            if (!std::isfinite(p_minus)) {
                return false;
            }
            dpdrho = (p_plus - p_minus) / (2.0 * h);
        }
        else {
            dpdrho = (p_plus - p_calc) / h;
        }
    }
    catch (const std::exception&) {
        return false;
    }
    if (!(std::isfinite(dpdrho) && dpdrho > 0.0)) {
        return false;
    }

    double gres = 0.0;
    try {
        gres = gres_cpp(t, rho, x, cppargs);
    }
    catch (const std::exception&) {
        return false;
    }
    if (!std::isfinite(gres)) {
        return false;
    }

    candidate->rho = rho;
    candidate->gres = gres;
    candidate->rel_resid = rel_resid;
    candidate->abs_p_error = abs_p_error;
    candidate->dpdrho = dpdrho;
    candidate->valid = true;
    return true;
}

bool density_root_from_seed_cpp(
    double t,
    double p,
    const vector<double> &x,
    int phase,
    const add_args &cppargs,
    double rho_seed,
    DensityRootCandidate *candidate,
    double *rho_root_out
) {
    if (!(std::isfinite(rho_seed) && rho_seed > 0.0)) {
        return false;
    }

    const int ncomp = static_cast<int>(x.size());
    const double nu_min = 1e-13;
    const double nu_max = 0.7405 - 1e-4;
    const double rho_min = reduced_density_to_molar(nu_min, t, ncomp, x, cppargs);
    const double rho_max = reduced_density_to_molar(nu_max, t, ncomp, x, cppargs);
    if (!(std::isfinite(rho_min) && std::isfinite(rho_max) && rho_max > rho_min)) {
        return false;
    }

    auto try_residual = [&](double rho, double *resid_out) -> bool {
        try {
            double resid = density_root_residual_cpp(rho, t, p, x, cppargs);
            if (!std::isfinite(resid)) {
                return false;
            }
            *resid_out = resid;
            return true;
        }
        catch (const std::exception&) {
            return false;
        }
    };

    rho_seed = std::max(rho_min, std::min(rho_seed, rho_max));
    double resid_seed = 0.0;
    if (try_residual(rho_seed, &resid_seed)) {
        DensityRootCandidate seed_candidate;
        if (density_root_valid_cpp(t, p, x, cppargs, rho_seed, &seed_candidate)) {
            if (candidate != nullptr) {
                *candidate = seed_candidate;
                candidate->rho_sort = rho_seed;
            }
            if (rho_root_out != nullptr) {
                *rho_root_out = rho_seed;
            }
            return true;
        }
    }

    double rho_lo = std::max(rho_min, rho_seed / 1.5);
    double rho_hi = std::min(rho_max, rho_seed * 1.5);
    double resid_lo = 0.0;
    double resid_hi = 0.0;
    bool have_lo = try_residual(rho_lo, &resid_lo);
    bool have_hi = try_residual(rho_hi, &resid_hi);

    const double expand_factor = 1.8;
    const int max_expansions = 14;
    bool bracketed = false;
    for (int iter = 0; iter < max_expansions; ++iter) {
        if (have_lo && have_hi && resid_lo * resid_hi < 0.0) {
            bracketed = true;
            break;
        }
        if (have_lo && try_residual(rho_seed, &resid_seed) && resid_lo * resid_seed < 0.0) {
            rho_hi = rho_seed;
            resid_hi = resid_seed;
            bracketed = true;
            break;
        }
        if (have_hi && try_residual(rho_seed, &resid_seed) && resid_seed * resid_hi < 0.0) {
            rho_lo = rho_seed;
            resid_lo = resid_seed;
            bracketed = true;
            break;
        }

        bool expanded = false;
        if (rho_lo > rho_min * (1.0 + 1e-12)) {
            double next_lo = std::max(rho_min, rho_lo / expand_factor);
            if (next_lo < rho_lo) {
                rho_lo = next_lo;
                have_lo = try_residual(rho_lo, &resid_lo);
                expanded = true;
            }
        }
        if (rho_hi < rho_max * (1.0 - 1e-12)) {
            double next_hi = std::min(rho_max, rho_hi * expand_factor);
            if (next_hi > rho_hi) {
                rho_hi = next_hi;
                have_hi = try_residual(rho_hi, &resid_hi);
                expanded = true;
            }
        }
        if (!expanded) {
            break;
        }
    }

    if (!bracketed) {
        return false;
    }

    try {
        double rho_root = density_brent_cpp(t, p, x, phase, cppargs, rho_lo, rho_hi, DBL_EPSILON, 1e-14, 200);
        DensityRootCandidate root_candidate;
        if (!density_root_valid_cpp(t, p, x, cppargs, rho_root, &root_candidate)) {
            return false;
        }
        root_candidate.rho_sort = rho_root;
        if (candidate != nullptr) {
            *candidate = root_candidate;
        }
        if (rho_root_out != nullptr) {
            *rho_root_out = rho_root;
        }
        return true;
    }
    catch (const std::exception&) {
        return false;
    }
}

DensityRootCandidate density_near_root_candidate_cpp(
    double t,
    double p,
    const vector<double> &x,
    const add_args &cppargs,
    double rho
) {
    DensityRootCandidate candidate;
    candidate.rho = rho;
    candidate.rho_sort = rho;
    candidate.abs_p_error = 1.0e300;
    candidate.rel_resid = 1.0e300;
    candidate.dpdrho = 1.0e300;
    candidate.gres = 1.0e300;
    try {
        double p_calc = p_cpp(t, rho, x, cppargs);
        if (std::isfinite(p_calc)) {
            candidate.abs_p_error = std::abs(p_calc - p);
            candidate.rel_resid = candidate.abs_p_error / std::max(std::abs(p), 1e-300);
        }
    }
    catch (const std::exception&) {
    }
    try {
        double h = std::max(1e-12, std::abs(rho) * 1e-6);
        double p_plus = p_cpp(t, rho + h, x, cppargs);
        if (std::isfinite(p_plus)) {
            if (rho > h) {
                double p_minus = p_cpp(t, rho - h, x, cppargs);
                if (std::isfinite(p_minus)) {
                    candidate.dpdrho = (p_plus - p_minus) / (2.0 * h);
                }
            }
            else {
                double p_calc = p_cpp(t, rho, x, cppargs);
                if (std::isfinite(p_calc)) {
                    candidate.dpdrho = (p_plus - p_calc) / h;
                }
            }
        }
    }
    catch (const std::exception&) {
    }
    try {
        double gres = gres_cpp(t, rho, x, cppargs);
        if (std::isfinite(gres)) {
            candidate.gres = gres;
        }
    }
    catch (const std::exception&) {
    }
    candidate.valid = std::isfinite(candidate.abs_p_error)
        && std::isfinite(candidate.dpdrho)
        && candidate.dpdrho > 0.0
        && std::isfinite(candidate.gres);
    return candidate;
}

DensityCandidateDiagnostics density_candidate_diagnostics_cpp(const DensityRootCandidate &candidate) {
    DensityCandidateDiagnostics out;
    out.rho_sort = candidate.rho_sort;
    out.rho = candidate.rho;
    out.gres = candidate.gres;
    out.rel_resid = candidate.rel_resid;
    out.abs_p_error = candidate.abs_p_error;
    out.dpdrho = candidate.dpdrho;
    out.valid = candidate.valid;
    return out;
}

DensitySolveResult density_solve_report_cpp(double t, double p, vector<double> x, int phase, const add_args &cppargs) {
    DensitySolveResult out;
    DensitySolveDiagnostics diagnostics;
    diagnostics.phase_kind = (phase == 1) ? "vap" : "liq";
    diagnostics.t = t;
    diagnostics.p = p;
    diagnostics.composition = x;
    diagnostics.validity_gate = "failed";
    diagnostics.best_near_root.abs_p_error = 1.0e300;
    diagnostics.best_near_root.rel_resid = 1.0e300;
    diagnostics.best_near_root.dpdrho = 1.0e300;
    diagnostics.best_near_root.gres = 1.0e300;

    int ncomp = static_cast<int>(x.size());
    vector<double> scan_grid = density_scan_grid_cpp();
    vector<DensityScanPoint> scan_points;
    scan_points.reserve(scan_grid.size());
    for (double nu : scan_grid) {
        DensityScanPoint point = density_scan_point_cpp(nu, t, ncomp, x, p, cppargs);
        scan_points.push_back(point);
        if (point.finite) {
            diagnostics.finite_point_count += 1;
            DensityRootCandidate near_candidate = density_near_root_candidate_cpp(t, p, x, cppargs, point.rho);
            if (near_candidate.abs_p_error < diagnostics.best_near_root.abs_p_error) {
                diagnostics.best_near_root = density_candidate_diagnostics_cpp(near_candidate);
            }
        }
    }
    diagnostics.scan_point_count = static_cast<int>(scan_points.size());

    vector<DensityBracket> coarse_brackets = density_brackets_cpp(scan_points);
    diagnostics.coarse_bracket_count = static_cast<int>(coarse_brackets.size());
    vector<DensityBracket> refined_brackets;
    for (const DensityBracket &coarse : coarse_brackets) {
        refine_density_brackets_cpp(coarse, t, ncomp, x, p, cppargs, refined_brackets);
    }
    diagnostics.refined_bracket_count = static_cast<int>(refined_brackets.size());

    if (refined_brackets.empty()) {
        diagnostics.fallback_used = true;
        if (std::isfinite(diagnostics.best_near_root.rho) && diagnostics.best_near_root.rho > 0.0) {
            DensityRootCandidate fallback_candidate;
            double rho_root = 0.0;
            if (density_root_from_seed_cpp(t, p, x, phase, cppargs, diagnostics.best_near_root.rho, &fallback_candidate, &rho_root)) {
                diagnostics.candidate_roots.push_back(density_candidate_diagnostics_cpp(fallback_candidate));
                diagnostics.candidate_root_count = static_cast<int>(diagnostics.candidate_roots.size());
                diagnostics.validity_gate = "passed";
                diagnostics.rejection_reason = "";
                out.rho = rho_root;
                out.valid = true;
                out.diagnostics = diagnostics;
                return out;
            }
        }
        diagnostics.fallback_rejected_reason = "no refined density brackets";
        diagnostics.rejection_reason = "No continuous density root brackets were found for the requested state";
        out.diagnostics = diagnostics;
        return out;
    }

    vector<DensityRootCandidate> candidates;
    candidates.reserve(refined_brackets.size());
    for (const DensityBracket &bracket : refined_brackets) {
        DensityRootCandidate candidate;
        candidate.rho_sort = ::reduced_density_to_molar(0.5 * (bracket.nu_lo + bracket.nu_hi), t, ncomp, x, cppargs);
        candidate.abs_p_error = 1.0e300;
        candidate.rel_resid = 1.0e300;
        candidate.dpdrho = 1.0e300;
        candidate.gres = 1.0e300;
        try {
            double rho_lo = ::reduced_density_to_molar(bracket.nu_lo, t, ncomp, x, cppargs);
            double rho_hi = ::reduced_density_to_molar(bracket.nu_hi, t, ncomp, x, cppargs);
            double rho_root = ::density_brent_cpp(t, p, x, phase, cppargs, rho_lo, rho_hi, DBL_EPSILON, 1e-14, 200);
            if (!density_root_valid_cpp(t, p, x, cppargs, rho_root, &candidate)) {
                candidate = density_near_root_candidate_cpp(t, p, x, cppargs, rho_root);
            }
            candidate.rho_sort = candidate.rho;
        }
        catch (const std::exception&) {
            candidate.valid = false;
        }
        if (candidate.abs_p_error < diagnostics.best_near_root.abs_p_error) {
            diagnostics.best_near_root = density_candidate_diagnostics_cpp(candidate);
        }
        candidates.push_back(candidate);
    }

    diagnostics.candidate_root_count = static_cast<int>(candidates.size());
    for (const DensityRootCandidate &candidate : candidates) {
        diagnostics.candidate_roots.push_back(density_candidate_diagnostics_cpp(candidate));
    }
    if (candidates.empty()) {
        diagnostics.rejection_reason = "Density solver did not produce any candidate roots";
        out.diagnostics = diagnostics;
        return out;
    }

    std::sort(candidates.begin(), candidates.end(), [](const DensityRootCandidate &a, const DensityRootCandidate &b) {
        return a.rho_sort < b.rho_sort;
    });

    const double rho_tol = 1e-8;
    DensityRootCandidate *best = nullptr;
    if (phase == 1) {
        const double rho_extreme = candidates.front().rho_sort;
        for (DensityRootCandidate &candidate : candidates) {
            if (std::abs(candidate.rho_sort - rho_extreme) > rho_tol * std::max(1.0, std::abs(rho_extreme))) {
                break;
            }
            if (candidate.valid && (best == nullptr || candidate.gres < best->gres)) {
                best = &candidate;
            }
        }
        diagnostics.rejection_reason = "No valid density root found for vapor phase";
    }
    else {
        const double rho_extreme = candidates.back().rho_sort;
        for (auto it = candidates.rbegin(); it != candidates.rend(); ++it) {
            if (std::abs(it->rho_sort - rho_extreme) > rho_tol * std::max(1.0, std::abs(rho_extreme))) {
                break;
            }
            if (it->valid && (best == nullptr || it->gres < best->gres)) {
                best = &(*it);
            }
        }
        diagnostics.rejection_reason = "No valid density root found for liquid phase";
    }
    if (best != nullptr) {
        out.rho = best->rho;
        out.valid = true;
        diagnostics.validity_gate = "passed";
        diagnostics.rejection_reason = "";
        out.diagnostics = diagnostics;
        return out;
    }
    diagnostics.fallback_used = true;
    if (std::isfinite(diagnostics.best_near_root.rho) && diagnostics.best_near_root.rho > 0.0) {
        DensityRootCandidate fallback_candidate;
        double rho_root = 0.0;
        if (density_root_from_seed_cpp(t, p, x, phase, cppargs, diagnostics.best_near_root.rho, &fallback_candidate, &rho_root)) {
            diagnostics.candidate_roots.push_back(density_candidate_diagnostics_cpp(fallback_candidate));
            diagnostics.candidate_root_count = static_cast<int>(diagnostics.candidate_roots.size());
            diagnostics.validity_gate = "passed";
            diagnostics.rejection_reason = "";
            out.rho = rho_root;
            out.valid = true;
            out.diagnostics = diagnostics;
            return out;
        }
    }
    diagnostics.fallback_rejected_reason = diagnostics.rejection_reason;
    out.diagnostics = diagnostics;
    return out;
}

double den_cpp(double t, double p, vector<double> x, int phase, const add_args &cppargs) {
    /**
    Solve for the molar density when temperature and pressure are given.

    Parameters
    ----------
    t : double
        Temperature (K)
    p : double
        Pressure (Pa)
    x : vector<double>, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    phase : int
        The phase for which the calculation is performed. Options: 0 (liquid),
        1 (vapor).
    cppargs : add_args
        A struct containing additional arguments that can be passed for
        use in PC-SAFT:

        m : vector<double>, shape (n,)
            Segment number for each component.
        s : vector<double>, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : vector<double>, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : vector<double>, shape (n*n,)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : vector<double>, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : vector<double>, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        z : vector<double>, shape (n,)
            Charge number of the ions
        dielc : double
            Dielectric constant of the medium to be used for electrolyte
            calculations.

    Returns
    -------
    rho : double
        Molar density (mol m^-3)
    */
    int ncomp = static_cast<int>(x.size());
    vector<double> scan_grid = density_scan_grid_cpp();
    vector<DensityScanPoint> scan_points;
    scan_points.reserve(scan_grid.size());
    for (double nu : scan_grid) {
        scan_points.push_back(density_scan_point_cpp(nu, t, ncomp, x, p, cppargs));
    }

    vector<DensityBracket> coarse_brackets = density_brackets_cpp(scan_points);
    vector<DensityBracket> refined_brackets;
    for (const DensityBracket &coarse : coarse_brackets) {
        refine_density_brackets_cpp(coarse, t, ncomp, x, p, cppargs, refined_brackets);
    }

    if (refined_brackets.empty()) {
        throw SolutionError(density_failure_message_cpp(
            "No continuous density root brackets were found for the requested state", t, p, x, phase));
    }

    vector<DensityRootCandidate> candidates;
    candidates.reserve(refined_brackets.size());
    for (const DensityBracket &bracket : refined_brackets) {
        DensityRootCandidate candidate;
        candidate.rho_sort = ::reduced_density_to_molar(0.5 * (bracket.nu_lo + bracket.nu_hi), t, ncomp, x, cppargs);

        try {
            double rho_lo = ::reduced_density_to_molar(bracket.nu_lo, t, ncomp, x, cppargs);
            double rho_hi = ::reduced_density_to_molar(bracket.nu_hi, t, ncomp, x, cppargs);
            double rho_root = ::density_brent_cpp(t, p, x, phase, cppargs, rho_lo, rho_hi, DBL_EPSILON, 1e-14, 200);
            density_root_valid_cpp(t, p, x, cppargs, rho_root, &candidate);
        }
        catch (const std::exception&) {
            candidate.valid = false;
        }

        candidates.push_back(candidate);
    }

    if (candidates.empty()) {
        throw SolutionError(density_failure_message_cpp(
            "Density solver did not produce any candidate roots", t, p, x, phase));
    }

    std::sort(candidates.begin(), candidates.end(), [](const DensityRootCandidate &a, const DensityRootCandidate &b) {
        return a.rho_sort < b.rho_sort;
    });

    const double rho_tol = 1e-8;
    if (phase == 1) {
        const double rho_extreme = candidates.front().rho_sort;
        DensityRootCandidate *best = nullptr;
        for (DensityRootCandidate &candidate : candidates) {
            if (std::abs(candidate.rho_sort - rho_extreme) > rho_tol * std::max(1.0, std::abs(rho_extreme))) {
                break;
            }
            if (candidate.valid && (best == nullptr || candidate.gres < best->gres)) {
                best = &candidate;
            }
        }
        if (best != nullptr) {
            return best->rho;
        }
        throw SolutionError(density_failure_message_cpp(
            "No valid density root found for vapor phase", t, p, x, phase));
    }

    const double rho_extreme = candidates.back().rho_sort;
    DensityRootCandidate *best = nullptr;
    for (auto it = candidates.rbegin(); it != candidates.rend(); ++it) {
        if (std::abs(it->rho_sort - rho_extreme) > rho_tol * std::max(1.0, std::abs(rho_extreme))) {
            break;
        }
        if (it->valid && (best == nullptr || it->gres < best->gres)) {
            best = &(*it);
        }
    }
    if (best != nullptr) {
        return best->rho;
    }
    throw SolutionError(density_failure_message_cpp(
        "No valid density root found for liquid phase", t, p, x, phase));
}

// EqID: rho_from_eta
// EqID: rho_reduced
double reduced_density_to_molar(double nu, double t, int ncomp, vector<double> x, const add_args &cppargs) {
    vector<double> d(ncomp);
    double summ = 0.;
    for (int i = 0; i < ncomp; i++) {
        d[i] = cppargs.s[i]*(1-0.12*std::exp(-3*cppargs.e[i] / t));
        if (!cppargs.z.empty() && is_ion_species(cppargs, i)) {
            d[i] = ion_diameter_cpp(i, t, cppargs);
        }
        summ += x[i]*cppargs.m[i]*pow(d[i],3.);
    }

    return 6/PI*nu/summ*1.0e30/N_AV;
}

/**
This function implements a 1-D bounded solver using the algorithm from Brent, R. P., Algorithms for Minimization Without Derivatives.
Englewood Cliffs, NJ: Prentice-Hall, 1973. Ch. 3-4.

a and b must bound the solution of interest and f(a) and f(b) must have opposite signs.  If the function is continuous, there must be
at least one solution in the interval [a,b].

@param a The minimum bound for the solution of f=0
@param b The maximum bound for the solution of f=0
@param macheps The machine precision
@param tol_abs Tolerance (absolute)
@param maxiter Maximum number of steps allowed.  Will throw a SolutionError if the solution cannot be found
*/
double density_brent_cpp(double t, double p, vector<double> x, int phase, const add_args &cppargs, double a, double b,
    double macheps, double tol_abs, int maxiter)
{
    int iter;
    double fa,fb,c,fc,m,tol,d,e,pp,q,s,r;
    fa = ::density_root_residual_cpp(a, t, p, x, cppargs);
    fb = ::density_root_residual_cpp(b, t, p, x, cppargs);

    // If one of the boundaries is to within tolerance, just stop
    if (std::abs(fb) < tol_abs) { return b;}
    if (std::isnan(fb)){
        throw ValueError("density root solver f(b) is NAN for b");
    }
    if (std::abs(fa) < tol_abs) { return a;}
    if (std::isnan(fa)){
        throw ValueError("density root solver f(a) is NAN for a");
    }
    if (fa*fb>0){
        throw ValueError("density root solver inputs do not bracket the root");
    }

    c=a;
    fc=fa;
    iter=1;
    if (std::abs(fc)<std::abs(fb)){
        // Goto ext: from Brent root solver ALGOL code
        a=b;
        b=c;
        c=a;
        fa=fb;
        fb=fc;
        fc=fa;
    }
    d=b-a;
    e=b-a;
    m=0.5*(c-b);
    tol=2*macheps*std::abs(b)+tol_abs;
    while (std::abs(m)>tol && fb!=0){
        // See if a bisection is forced
        if (std::abs(e)<tol || std::abs(fa) <= std::abs(fb)){
            m=0.5*(c-b);
            d=e=m;
        }
        else{
            s=fb/fa;
            if (a==c){
                //Linear interpolation
                pp=2*m*s;
                q=1-s;
            }
            else{
                //Inverse quadratic interpolation
                q=fa/fc;
                r=fb/fc;
                m=0.5*(c-b);
                pp=s*(2*m*q*(q-r)-(b-a)*(r-1));
                q=(q-1)*(r-1)*(s-1);
            }
            if (pp>0){
                q=-q;
            }
            else{
                pp=-pp;
            }
            s=e;
            e=d;
            m=0.5*(c-b);
            if (2*pp<3*m*q-std::abs(tol*q) || pp<std::abs(0.5*s*q)){
                d=pp/q;
            }
            else{
                m=0.5*(c-b);
                d=e=m;
            }
        }
        a=b;
        fa=fb;
        if (std::abs(d)>tol){
            b+=d;
        }
        else if (m>0){
            b+=tol;
        }
        else{
            b+=-tol;
        }
        fb=::density_root_residual_cpp(b, t, p, x, cppargs);
        if (std::isnan(fb)){
            throw ValueError("density root solver f(t) is NAN for t");
        }
        if (std::abs(fb) < macheps){
            return b;
        }
        if (fb*fc>0){
            c=a;
            fc=fa;
            d=e=b-a;
        }
        if (std::abs(fc)<std::abs(fb)){
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        m=0.5*(c-b);
        tol=2*macheps*std::abs(b)+tol_abs;
        iter+=1;
        if (std::isnan(a)){
            throw ValueError("density root solver a is NAN");}
        if (std::isnan(b)){
            throw ValueError("density root solver b is NAN");}
        if (std::isnan(c)){
            throw ValueError("density root solver c is NAN");}
        if (iter>maxiter){
            throw SolutionError("density root solver reached maximum number of steps");}
        if (std::abs(fb)< 2*macheps*std::abs(b)){
            return b;
        }
    }
    return b;
}

// EqID: density_solve_residual
double density_root_residual_cpp(double rhomolar, double t, double p, vector<double> x, const add_args &cppargs){
    double peos = p_cpp(t, rhomolar, x, cppargs);
    double pressure_scale = std::max(std::abs(p), 1e-3);
    double cost = (peos-p)/pressure_scale;
    if (std::isfinite(cost)) {
        return cost;
    }
    else {
        return HUGE_DBL;
    }
}
