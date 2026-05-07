# Graph Report - .  (2026-05-07)

## Corpus Check
- cluster-only mode — file stats not available

## Summary
- 1965 nodes · 4809 edges · 106 communities (91 shown, 15 thin omitted)
- Extraction: 85% EXTRACTED · 15% INFERRED · 0% AMBIGUOUS · INFERRED: 740 edges (avg confidence: 0.77)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `378edd7f`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- [[_COMMUNITY_Community 0|Community 0]]
- [[_COMMUNITY_Community 1|Community 1]]
- [[_COMMUNITY_Community 2|Community 2]]
- [[_COMMUNITY_Community 3|Community 3]]
- [[_COMMUNITY_Community 4|Community 4]]
- [[_COMMUNITY_Community 5|Community 5]]
- [[_COMMUNITY_Community 6|Community 6]]
- [[_COMMUNITY_Community 7|Community 7]]
- [[_COMMUNITY_Community 8|Community 8]]
- [[_COMMUNITY_Community 9|Community 9]]
- [[_COMMUNITY_Community 10|Community 10]]
- [[_COMMUNITY_Community 11|Community 11]]
- [[_COMMUNITY_Community 12|Community 12]]
- [[_COMMUNITY_Community 13|Community 13]]
- [[_COMMUNITY_Community 14|Community 14]]
- [[_COMMUNITY_Community 15|Community 15]]
- [[_COMMUNITY_Community 16|Community 16]]
- [[_COMMUNITY_Community 17|Community 17]]
- [[_COMMUNITY_Community 18|Community 18]]
- [[_COMMUNITY_Community 19|Community 19]]
- [[_COMMUNITY_Community 20|Community 20]]
- [[_COMMUNITY_Community 21|Community 21]]
- [[_COMMUNITY_Community 22|Community 22]]
- [[_COMMUNITY_Community 23|Community 23]]
- [[_COMMUNITY_Community 24|Community 24]]
- [[_COMMUNITY_Community 25|Community 25]]
- [[_COMMUNITY_Community 26|Community 26]]
- [[_COMMUNITY_Community 27|Community 27]]
- [[_COMMUNITY_Community 28|Community 28]]
- [[_COMMUNITY_Community 29|Community 29]]
- [[_COMMUNITY_Community 30|Community 30]]
- [[_COMMUNITY_Community 31|Community 31]]
- [[_COMMUNITY_Community 33|Community 33]]
- [[_COMMUNITY_Community 34|Community 34]]
- [[_COMMUNITY_Community 35|Community 35]]
- [[_COMMUNITY_Community 36|Community 36]]
- [[_COMMUNITY_Community 37|Community 37]]
- [[_COMMUNITY_Community 38|Community 38]]
- [[_COMMUNITY_Community 39|Community 39]]
- [[_COMMUNITY_Community 40|Community 40]]
- [[_COMMUNITY_Community 41|Community 41]]
- [[_COMMUNITY_Community 42|Community 42]]
- [[_COMMUNITY_Community 43|Community 43]]
- [[_COMMUNITY_Community 44|Community 44]]
- [[_COMMUNITY_Community 45|Community 45]]
- [[_COMMUNITY_Community 46|Community 46]]
- [[_COMMUNITY_Community 47|Community 47]]
- [[_COMMUNITY_Community 48|Community 48]]
- [[_COMMUNITY_Community 49|Community 49]]
- [[_COMMUNITY_Community 50|Community 50]]
- [[_COMMUNITY_Community 51|Community 51]]
- [[_COMMUNITY_Community 52|Community 52]]
- [[_COMMUNITY_Community 53|Community 53]]
- [[_COMMUNITY_Community 54|Community 54]]
- [[_COMMUNITY_Community 55|Community 55]]
- [[_COMMUNITY_Community 56|Community 56]]
- [[_COMMUNITY_Community 57|Community 57]]
- [[_COMMUNITY_Community 59|Community 59]]
- [[_COMMUNITY_Community 60|Community 60]]
- [[_COMMUNITY_Community 61|Community 61]]
- [[_COMMUNITY_Community 62|Community 62]]
- [[_COMMUNITY_Community 63|Community 63]]
- [[_COMMUNITY_Community 64|Community 64]]
- [[_COMMUNITY_Community 65|Community 65]]
- [[_COMMUNITY_Community 66|Community 66]]
- [[_COMMUNITY_Community 67|Community 67]]
- [[_COMMUNITY_Community 68|Community 68]]
- [[_COMMUNITY_Community 69|Community 69]]
- [[_COMMUNITY_Community 70|Community 70]]
- [[_COMMUNITY_Community 71|Community 71]]
- [[_COMMUNITY_Community 72|Community 72]]
- [[_COMMUNITY_Community 73|Community 73]]
- [[_COMMUNITY_Community 74|Community 74]]
- [[_COMMUNITY_Community 75|Community 75]]
- [[_COMMUNITY_Community 76|Community 76]]
- [[_COMMUNITY_Community 77|Community 77]]
- [[_COMMUNITY_Community 78|Community 78]]
- [[_COMMUNITY_Community 79|Community 79]]
- [[_COMMUNITY_Community 82|Community 82]]
- [[_COMMUNITY_Community 93|Community 93]]
- [[_COMMUNITY_Community 94|Community 94]]
- [[_COMMUNITY_Community 95|Community 95]]
- [[_COMMUNITY_Community 96|Community 96]]
- [[_COMMUNITY_Community 97|Community 97]]
- [[_COMMUNITY_Community 98|Community 98]]
- [[_COMMUNITY_Community 99|Community 99]]
- [[_COMMUNITY_Community 100|Community 100]]
- [[_COMMUNITY_Community 101|Community 101]]
- [[_COMMUNITY_Community 102|Community 102]]
- [[_COMMUNITY_Community 103|Community 103]]
- [[_COMMUNITY_Community 104|Community 104]]
- [[_COMMUNITY_Community 105|Community 105]]

## God Nodes (most connected - your core abstractions)
1. `main()` - 186 edges
2. `ValueError()` - 159 edges
3. `InputError` - 126 edges
4. `plot_figure_1()` - 70 edges
5. `ePCSAFTState` - 56 edges
6. `ePCSAFTMixture` - 38 edges
7. `epcsaft_density()` - 37 edges
8. `_fit_pure_ion_internal()` - 30 edges
9. `_fit_binary_pair_internal()` - 30 edges
10. `get_prop_dict()` - 28 edges

## Surprising Connections (you probably didn't know these)
- `evaluate_generic_regression_derivatives()` --calls--> `test_public_generic_derivative_evaluator_exposes_jacobian_and_hessian_skeleton()`  [INFERRED]
  src/epcsaft/regression.py → tests/api/test_regression_api.py
- `run_confidence_suite()` --calls--> `test_confidence_suite_smoke_mode_writes_bounded_report()`  [INFERRED]
  src/epcsaft/equilibrium_core/confidence.py → tests/equilibrium/test_electrolyte_lle_confidence.py
- `run_confidence_suite()` --calls--> `test_opt_in_confidence_report_generates_full_outputs()`  [INFERRED]
  src/epcsaft/equilibrium_core/confidence.py → tests/equilibrium/test_electrolyte_lle_confidence.py
- `run_script()` --calls--> `main()`  [EXTRACTED]
  analyses/2012_held/scripts/run_all.py → src/epcsaft/__main__.py
- `save_figure()` --calls--> `paper_validation_output_path()`  [INFERRED]
  analyses/2012_held/scripts/_common.py → scripts/plot_outputs.py

## Communities (106 total, 15 thin omitted)

### Community 0 - "Community 0"
Cohesion: 0.05
Nodes (94): Return the osmotic coefficient., _fit_pure_neutral_least_squares_internal(), Internal native least-squares comparison hook for pure-neutral regression., activity_coefficient_values_cpp(), assign_activity_aux_cpp(), assign_activity_metadata_cpp(), assign_pair_activity_cpp(), build_infinite_dilution_reference_cpp() (+86 more)

### Community 1 - "Community 1"
Cohesion: 0.05
Nodes (91): association_intermediate_state_cpp(), association_setup_cpp(), association_site_fraction_composition_terms_cpp(), association_site_fractions_cpp(), dadrho_assoc_cpp(), dadt_assoc_cpp(), dadx_assoc_cpp(), solve_association_site_fractions_cpp() (+83 more)

### Community 2 - "Community 2"
Cohesion: 0.05
Nodes (58): build_workbook(), _infinite_dilution_state(), _write_active_mode_detail(), _write_active_mode_summary(), _effective_diameters(), _state_for_ion(), _calc_ln_miac_contributions(), _inf_dilution_state() (+50 more)

### Community 3 - "Community 3"
Cohesion: 0.05
Nodes (57): _minimal_nacl_records(), _stub_native_generic_runner(), test_fit_binary_pair_can_fit_all_constant_binary_interaction_targets(), test_fit_binary_pair_vle_kij_default_and_rejects_temperature_models(), test_fit_pure_ion_accepts_d_born_and_born_user_options(), test_fit_pure_ion_default_s_e_bounds_and_multistart_contract(), test_fit_pure_ion_passes_explicit_mean_ionic_pair_label_to_native_backend(), test_fit_pure_neutral_rejects_non_phase1_targets() (+49 more)

### Community 4 - "Community 4"
Cohesion: 0.08
Nodes (70): solve_equilibrium_native_binding(), _phase_state(), clip_normalize(), component_rich_composition(), composition_charge(), composition_from_log_weights(), damping_schedule(), l2_norm() (+62 more)

### Community 5 - "Community 5"
Cohesion: 0.07
Nodes (52): array_to_double_vector(), array_to_int_vector(), chemical_options_from_request(), density_records_from_arrays(), electrolyte_bubble_options_from_request(), generic_record_from_dict(), generic_records_from_list(), generic_regression_debug_to_dict() (+44 more)

### Community 6 - "Community 6"
Cohesion: 0.07
Nodes (51): _aad_summary(), add_figure_caption(), baygi_output_root(), _baygi_regression_bounds(), build_feed_formula_from_salt_free_molefractions(), build_feed_formula_from_total_feed_weights(), closest_group_to_weight_fraction(), configure_style() (+43 more)

### Community 7 - "Community 7"
Cohesion: 0.08
Nodes (57): EquilibriumOptions, Numerical controls for equilibrium solvers., _write_csv(), _normalize_pair(), Exception raised when a solver fails to converge or returns invalid data., SolutionError, _all_tielines_plot(), _bar_comparison_plot() (+49 more)

### Community 8 - "Community 8"
Cohesion: 0.08
Nodes (53): _charge(), compare_khudaida_aad_tables(), compare_khudaida_digitized_paper_to_package(), evaluate_khudaida_solver_gate(), evaluate_khudaida_tieline(), _feed_row_to_formula(), _figure_data_dir(), _finite_or_none() (+45 more)

### Community 9 - "Community 9"
Cohesion: 0.07
Nodes (38): _figure3_paper_values(), _summary(), _curve_bundle(), _interp_metrics(), _read_digitized_series(), _curve_with_override(), _panel_map(), _plot_scan() (+30 more)

### Community 10 - "Community 10"
Cohesion: 0.1
Nodes (42): _as_rule_number(), available_datasets(), _canonical_solvent_tokens(), _coerce_bool(), _convert_weight_to_mole_fractions(), _deep_update(), default_user_options(), _deterministic_default() (+34 more)

### Community 11 - "Community 11"
Cohesion: 0.07
Nodes (40): _auto_ylim(), _composition_from_x_il(), _compute_gsolv(), _computed_row(), _curve_for_mode(), _figure2_params(), _ion_diameter_angstrom(), _load_optional_digitized() (+32 more)

### Community 12 - "Community 12"
Cohesion: 0.12
Nodes (37): _calc_curve(), comp_signature(), mw_mix(), normalized_comp(), salt_mole_fraction_from_molality(), stoich_for_salt(), _available_solvent_systems(), build_params_for_variant() (+29 more)

### Community 13 - "Community 13"
Cohesion: 0.11
Nodes (33): Solve a neutral liquid-liquid TP flash., Run transformed-basis electrolyte TP stability analysis., Run a native-backed equilibrium calculation for this mixture., Solve an ordered equilibrium curve, reusing each accepted split as the next seed, _add_legacy_option_diagnostics(), _call_native_equilibrium(), _diagnostics_with_legacy_candidate(), electrolyte_feed_from_molality() (+25 more)

### Community 14 - "Community 14"
Cohesion: 0.09
Nodes (18): test_explicit_finite_difference_dadx_reports_finite_difference_backend(), test_from_params_rejects_legacy_electrolyte_keys(), test_hc_dadx_autodiff_matches_analytic_terms(), test_ionic_activity_and_solution_methods_return_expected_values(), test_mu_born_autodiff_matches_analytic_born_terms(), test_mu_dh_autodiff_matches_analytic_ion_terms(), test_neutral_composition_and_fugacity_terms_return_expected_values(), test_neutral_scalar_methods_return_expected_values() (+10 more)

### Community 15 - "Community 15"
Cohesion: 0.12
Nodes (29): output_path(), _curve_for_combo(), _load_rows_for_combo(), make_plot_for_salt(), _single_solvent_combo_map(), analysis_final_dir(), analysis_final_path(), _analysis_root_for() (+21 more)

### Community 16 - "Community 16"
Cohesion: 0.08
Nodes (24): ePCSAFTState, ncomp(), Resolve supported short scientific aliases for state methods., Clear internal runtime caches used for repeated state/reference evaluations., Reset runtime cache hit/fallback counters without clearing cached values., Return a short debugging representation of the mixture., Immutable thermodynamic state bound to one mixture., Return the canonical full-name -> abbreviation map for state methods. (+16 more)

### Community 17 - "Community 17"
Cohesion: 0.12
Nodes (25): _apply_weights(), _best_scalar_weight(), _load_bookkeeping(), Scan alternative Figure 6b accounting methods against digitized contribution lin, _rmse(), run_scan(), _interp(), Generate Figure 6b fit-check plots from bookkeeping curves plus z-correction wei (+17 more)

### Community 18 - "Community 18"
Cohesion: 0.1
Nodes (26): test_reactive_electrolyte_bubble_accepts_phase_handoff_speciation_residuals(), test_reactive_electrolyte_bubble_respects_configured_phase_handoff_tolerances(), test_reactive_electrolyte_bubble_sweep_preserves_phase_handoff_tolerances(), electrolyte_bubble_pressure(), ElectrolyteBubbleResult, _normalize_species_labels(), _phase_by_label(), Native-only electrolyte bubble-pressure public contracts. (+18 more)

### Community 19 - "Community 19"
Cohesion: 0.16
Nodes (28): _assoc_is_enabled(), _build_single_component_params(), _coerce_bounds(), _ensure_native_vector_payload(), _family_scale(), _fit_derivative_metadata(), _fit_mea_co2_h2o_component(), _fit_pure_ion_internal() (+20 more)

### Community 20 - "Community 20"
Cohesion: 0.12
Nodes (18): test_solve_reactive_staged_equilibrium_returns_chemical_and_phase_results(), _assert_json_like(), _hydrocarbon_mixture(), _ionic_mixture(), test_equilibrium_phase_exposes_explicit_ln_fugacity_alias(), test_equilibrium_rejects_invalid_public_inputs(), test_equilibrium_rejects_ionic_mixtures_for_v1(), test_explicit_flash_tp_matches_legacy_equilibrium_dispatch() (+10 more)

### Community 21 - "Community 21"
Cohesion: 0.13
Nodes (22): ePCSAFTMixture, Return native runtime cache counters for profiling and validation., Solve a neutral TP flash with explicit thermodynamic-method naming., Solve a neutral bubble pressure from liquid composition at fixed temperature., Solve a neutral dew pressure from vapor composition at fixed temperature., Run neutral TP tangent-plane-distance stability analysis., Solve charge-constrained electrolyte LLE at fixed temperature and pressure., Solve an ordered equilibrium sweep with workflow-specific continuation. (+14 more)

### Community 22 - "Community 22"
Cohesion: 0.15
Nodes (21): _ascani_water_butanol_nacl_mixture(), _case2_feed(), _case2_mixture(), test_confidence_suite_smoke_mode_writes_bounded_report(), test_opt_in_confidence_report_generates_full_outputs(), test_ascani_case2_mixed_salt_solves_without_local_model_fixture(), test_ascani_counterion_basis_has_expected_rank_and_preserves_charge(), test_auto_kind_routes_explicit_ionic_feed_to_electrolyte_lle() (+13 more)

### Community 23 - "Community 23"
Cohesion: 0.11
Nodes (23): _benchmark_vector_map(), _best_pair_label(), _choose_pure_file(), evaluate_generic_regression_derivatives(), _find_matching_pure_files(), fit_pure_ion(), FitParameter, _native_target_payload() (+15 more)

### Community 24 - "Community 24"
Cohesion: 0.12
Nodes (22): Return the compressibility factor., Return the residual Helmholtz energy., Return the residual enthalpy., Return the residual entropy., Return the residual Gibbs energy., Return a diagnostic dictionary of the main state properties., Return the state density in the requested units., Return the molar density of the bound state in mol/m^3. (+14 more)

### Community 25 - "Community 25"
Cohesion: 0.15
Nodes (20): assert_plot_with_data(), save_comparison_plot(), test_electrolyte_lle_confidence_plots_are_written_to_local_output(), _ascani_electrolyte_lle_result(), test_equilibrium_electrolyte_khudaida_seeded_solver_gate_plot(), test_equilibrium_electrolyte_lle_phase_composition_plot(), test_equilibrium_electrolyte_lle_residual_closure_plot(), test_equilibrium_lle_tie_line_plot_is_written_to_gallery() (+12 more)

### Community 26 - "Community 26"
Cohesion: 0.13
Nodes (20): Solve homogeneous activity-coupled chemical equilibrium., Run staged chemical equilibrium followed by an explicit phase route., _finite_float_option(), _json_like(), _normalize_options(), _named_reaction_residuals(), _normalize_balances(), _normalize_composition() (+12 more)

### Community 27 - "Community 27"
Cohesion: 0.19
Nodes (16): _fit_line(), _model_state(), _safe_corr(), _assoc_breakdown(), _disp_breakdown(), _generic_contribution_rows(), _hc_breakdown(), _kjmol() (+8 more)

### Community 28 - "Community 28"
Cohesion: 0.18
Nodes (21): _binary_species_from_records(), _build_binary_terms(), _build_mea_co2_h2o_terms(), _build_pure_ion_terms(), _build_pure_neutral_terms(), _composition_from_record(), _infer_species_union(), _ion_composition_from_record() (+13 more)

### Community 29 - "Community 29"
Cohesion: 0.18
Nodes (19): _closest_group(), compare_all(), _compare_values(), _composition_json(), _figure9_curve_m_max(), _format_value(), generate_figure_4(), _group_by_weight() (+11 more)

### Community 30 - "Community 30"
Cohesion: 0.18
Nodes (20): _canonical_salt(), _composition_columns(), _format(), _lookup(), merge_and_sync(), _parse_salt_from_stem(), Merge/sync MIAC CSV variants (molality, miac_m, mole_fraction, miac).  This sc, _read_rows() (+12 more)

### Community 31 - "Community 31"
Cohesion: 0.26
Nodes (20): _assert_runtime_sentinels(), _clean_cell(), _empty_matrix(), _enforce_symmetry(), _extract_2005(), _extract_2014(), _extract_2020(), _extract_for_paper() (+12 more)

### Community 33 - "Community 33"
Cohesion: 0.15
Nodes (18): attach_code_refs(), check_matches(), docs_only_entries(), enforce_traceability(), is_documentation_only(), missing_cpp_ref_entries(), parse_code_refs(), parse_equations() (+10 more)

### Community 34 - "Community 34"
Cohesion: 0.17
Nodes (18): _baygi_mea_records(), _perturbed(), test_baygi_table2_mea_association_scheme_optimizer_does_not_worsen_seed(), test_baygi_table2_mea_parameters_score_better_than_perturbed_seed(), baygi_mea_density_molar(), baygi_mea_fit_records(), baygi_mea_psat_pa(), _bisect_log_pressure_root() (+10 more)

### Community 35 - "Community 35"
Cohesion: 0.16
Nodes (14): _clean_dist(), _env(), _newest_wheel(), _smoke_wheel(), _capture(), _clean(), _configure(), _configured_generator() (+6 more)

### Community 36 - "Community 36"
Cohesion: 0.22
Nodes (15): append_payload_rows(), _comparison_size(), _finish_figure(), math_label(), _max_wrapped_lines(), _maybe_use_symmetric_log_scale(), save_contribution_closure_plot(), save_contribution_term_breakdown_plot() (+7 more)

### Community 37 - "Community 37"
Cohesion: 0.14
Nodes (17): EquilibriumPhase, EquilibriumResult, _native_phase_from_payload(), _native_result_from_payload(), _native_trial_from_payload(), _neutral_bubble_dew_result(), Structured result returned by an equilibrium calculation., Return a JSON-like result payload. (+9 more)

### Community 38 - "Community 38"
Cohesion: 0.2
Nodes (13): _payload(), _legend_handles(), _plot_binary_reference(), _save_combined(), _weight_to_mole_comp(), _plot_panel(), _plot_40wt_solvent(), _series_xy() (+5 more)

### Community 39 - "Community 39"
Cohesion: 0.21
Nodes (17): _miac_curve(), build_params(), _freeze_comp(), _freeze_user_options(), _gsolv_ion_cached(), mean_ionic_activity_curve(), molality_to_species_molefraction(), pair_key_for_salt() (+9 more)

### Community 40 - "Community 40"
Cohesion: 0.18
Nodes (10): _log_k_from_state(), _mea_like_mixture(), _neutral_log_k_from_fugacity_activity(), test_native_chemical_equilibrium_matches_activity_coupled_salt_speciation(), test_native_chemical_equilibrium_soft_start_reports_for_hard_ionic_speciation(), test_native_chemical_equilibrium_solution_shifts_when_fugacity_model_changes(), test_native_chemical_equilibrium_solves_hard_mea_like_speciation_and_returns_phase_handoff(), test_native_chemical_equilibrium_uses_convex_soft_start_for_bad_neutral_seed() (+2 more)

### Community 41 - "Community 41"
Cohesion: 0.2
Nodes (16): _explicit_beta_to_formula_beta(), _explicit_to_formula_composition(), _formula_to_explicit_composition(), _ion_stem(), _independent_counterion_pairs(), _pair_for_salt_label(), _pair_payload(), Ascani-style transformed variables for electrolyte liquid splits. (+8 more)

### Community 42 - "Community 42"
Cohesion: 0.21
Nodes (16): _build_catalog(), _eval_csv_rows(), _linear_t(), Build ePC-SAFT electrolyte parameter catalog JSON + CSV reference tables.  Thi, _water_sigma_expr(), _dataset_binary_dir(), _dataset_dir(), _dataset_pure_dir() (+8 more)

### Community 43 - "Community 43"
Cohesion: 0.2
Nodes (14): classify_issue(), _clean_field_value(), fetch_issue(), _has_value(), issue_number_from_token(), IssueTriage, _labels_from_issue(), parse_issue_fields() (+6 more)

### Community 44 - "Community 44"
Cohesion: 0.21
Nodes (14): _build_strategy2_params(), _base_params(), calc_osmotic_curve(), _elec_model_no_born_constant(), load_osmotic_data(), mole_fraction_from_molality_11(), osmotic_molality_from_fugacity(), water_sigma() (+6 more)

### Community 45 - "Community 45"
Cohesion: 0.29
Nodes (14): initial_phases_from_result(), Build electrolyte LLE ``initial_phases`` from an accepted aq/org result., test_equilibrium_curve_uses_previous_hubach_split_as_seed(), test_initial_phases_from_result_round_trips_hubach_split(), _hubach_mixture(), _row0_feed(), _row0_initial_phases(), test_hubach_cold_start_density_failure_payload_is_json_safe() (+6 more)

### Community 46 - "Community 46"
Cohesion: 0.23
Nodes (16): from_dataset(), _apply_constant_mixed_rel_perm_precompute(), _apply_mixed_solvent_ion_dispersion(), _apply_mixed_solvent_ion_sigma(), _as_composition_array(), _compute_constant_mixed_rel_perm(), _compute_constant_salt_free_weight_avg_rel_perm(), get_prop_dict() (+8 more)

### Community 47 - "Community 47"
Cohesion: 0.19
Nodes (14): test_create_parameter_template_creates_loadable_scaffold(), test_write_fit_result_updates_ion_row_and_binary_matrix_symmetrically(), test_write_fit_result_updates_only_target_pure_row(), create_parameter_template(), _infer_pure_template_name(), _prompt(), Helpers for creating user-owned ePC-SAFT parameter templates., Create a user-owned dataset scaffold and return its root path.      If any of (+6 more)

### Community 48 - "Community 48"
Cohesion: 0.26
Nodes (14): Return the residual chemical potentials., _vector_terms_dict(), Convert a vector-like payload into a NumPy float array., vector_to_array(), mu_assoc_cpp(), mu_born_cpp(), mu_contribution_cpp(), mu_disp_cpp() (+6 more)

### Community 49 - "Community 49"
Cohesion: 0.23
Nodes (12): capabilities(), _direct_url_payload(), _git_commit(), _mtime_utc(), _native_extension_path(), _path_from_file_url(), Runtime metadata and capability discovery for downstream applications., Return JSON-like package, source, and native-extension metadata. (+4 more)

### Community 50 - "Community 50"
Cohesion: 0.22
Nodes (10): load_experimental_rows(), _objective(), _solve_lle(), _split_nh4cl_mass(), _candidate_formula_feeds(), formula_to_ion_basis(), ion_to_formula_basis(), _model_objective() (+2 more)

### Community 51 - "Community 51"
Cohesion: 0.21
Nodes (7): annotate_percent_deltas(), load_indexed_csv(), percent_delta(), _round_down_to_multiple(), _round_up_to_multiple(), set_strict_bar_ylim(), Table

### Community 52 - "Community 52"
Cohesion: 0.18
Nodes (9): _electrolyte_initial_phase_seed(), Public package interface for the native ePC-SAFT runtime., classify_equilibrium_route(), Route classification for Python equilibrium calculations., Return the internal route for a public equilibrium request., charge_neutral_lle_seed_from_org_phase(), Helpers for constructing charge-neutral electrolyte LLE initial phases., solvent_endpoint_seed() (+1 more)

### Community 53 - "Community 53"
Cohesion: 0.15
Nodes (13): _associating_pure_payload(), _benchmark_seed_payloads(), _debug_native_pure_neutral_objective(), evaluate_pure_neutral_derivatives(), _fit_mea_co2_h2o_pure_parameter_benchmark(), load_regression_records(), _normalize_fit_targets(), _normalize_records() (+5 more)

### Community 54 - "Community 54"
Cohesion: 0.15
Nodes (11): Build the full activity-coefficient payload for internal reuse., Return activity coefficients in the requested form., Return ion solvation free-energy values keyed by species., _resolve_solvent_override(), ActivityCoefficientResult, Return component activity coefficients keyed by species label., Return ion solvation free-energy values keyed by ion label., Return mean-ionic activity coefficients in the requested basis. (+3 more)

### Community 55 - "Community 55"
Cohesion: 0.33
Nodes (12): Return fugacity coefficients, defaulting to natural-log form., fugacity_coefficient(), lnfug_assoc_cpp(), lnfug_born_cpp(), lnfug_contribution_cpp(), lnfug_correction_scale_cpp(), lnfug_cpp(), lnfug_disp_cpp() (+4 more)

### Community 56 - "Community 56"
Cohesion: 0.39
Nodes (10): _methanol_cyclohexane_lle_benchmark(), test_lle_flash_distinct_stalled_seed_raises_solution_error(), test_lle_flash_phase_diagnostics_are_json_serializable_when_requested(), test_lle_flash_rejects_invalid_options_through_public_api(), test_lle_flash_rejects_invalid_public_inputs(), test_lle_flash_reports_no_split_for_identical_initial_phases(), test_lle_flash_without_initial_phases_finds_methanol_cyclohexane_split(), test_methanol_cyclohexane_lle_flash_closes_material_and_fugacity_balance() (+2 more)

### Community 57 - "Community 57"
Cohesion: 0.18
Nodes (12): _binary_seed_value(), _charge_map(), _copy_mapping(), fit_binary_pair(), _fit_binary_pair_internal(), _normalize_temperature_model(), _params_for_native_record(), Fit V1 binary interaction parameters against VLE x/y records. (+4 more)

### Community 59 - "Community 59"
Cohesion: 0.25
Nodes (6): _electrolyte_mixture(), test_native_electrolyte_stability_entrypoint_runs_in_cpp(), test_native_electrolyte_stability_honors_explicit_max_iterations(), test_public_electrolyte_lle_uses_native_backend_with_initial_phases(), test_public_electrolyte_stability_uses_native_backend(), test_public_equilibrium_result_comes_from_native_backend()

### Community 60 - "Community 60"
Cohesion: 0.27
Nodes (9): build_wheel(), _external_temp_root(), _has_build_dir(), _is_under(), _isolated_build_config(), PEP 517 backend wrapper for sandbox-safe Windows package builds., _sandbox_safe_mkdtemp(), _source_root() (+1 more)

### Community 61 - "Community 61"
Cohesion: 0.25
Nodes (10): test_runtime_options_accept_autodiff_modes_and_preserve_explicit_overrides(), test_runtime_options_default_to_auto_derivative_policy(), test_runtime_options_reject_legacy_electrolyte_shorthand(), test_runtime_options_reject_removed_polar_model(), _build_user_options_payload(), minimize_user_options(), _prune_default_overrides(), Normalize user options into the canonical runtime model schema. (+2 more)

### Community 62 - "Community 62"
Cohesion: 0.31
Nodes (6): _salt_speciation_mixture(), test_solve_reactive_speciation_best_effort_returns_nonconverged_result(), test_solve_reactive_speciation_concentration_standard_state_uses_molar_density(), test_solve_reactive_speciation_rejects_invalid_chemistry_inputs(), test_solve_reactive_speciation_returns_balanced_activity_coupled_state(), test_solve_reactive_speciation_strict_failure_reports_best_state()

### Community 63 - "Community 63"
Cohesion: 0.39
Nodes (8): _aard_by_series(), _artist_labels(), _diagnostic_rows(), _generate_baygi_outputs(), _seed_cached_baygi_outputs(), test_2015_baygi_figure_2_and_3_workflow_outputs_expected_series(), test_2015_baygi_figures_report_numeric_fit_quality(), metric_rows()

### Community 64 - "Community 64"
Cohesion: 0.31
Nodes (7): save_parity_plot(), test_runtime_alias_canonical_method_parity_plot(), test_runtime_diagnostics_public_method_parity_plot(), test_native_branch_and_contribution_reference_comparison_plot(), test_native_composition_derivative_finite_difference_parity_plot(), test_native_temperature_derivative_finite_difference_parity_plot(), test_runtime_pressure_density_constructor_parity_plot()

### Community 65 - "Community 65"
Cohesion: 0.5
Nodes (8): fit_mea_co2_h2o_electrolyte(), Fit the dataset-driven MEA-CO2-H2O electrolyte pure-parameter benchmark., _benchmark_dataset(), _benchmark_records(), test_mea_co2_h2o_benchmark_opt_in_real_multistart(), test_mea_co2_h2o_benchmark_smoke_fits_only_pure_parameters(), test_mea_co2_h2o_benchmark_writes_only_pure_rows_and_protects_existing_cells(), _write_pure_rows()

### Community 66 - "Community 66"
Cohesion: 0.22
Nodes (8): One reaction residual definition for reactive speciation., ReactionDefinition, _binary_interaction_declarations(), BinaryInteraction, FitTerm, One binary-interaction parameter requested for provenance-aware fitting., One weighted family of regression residuals., _source_token()

### Community 67 - "Community 67"
Cohesion: 0.33
Nodes (4): _body(), _issue_payload(), test_complete_downstream_issue_returns_ready_triage(), test_missing_minimal_reproducer_returns_needs_repro_recommendation()

### Community 68 - "Community 68"
Cohesion: 0.46
Nodes (7): _read(), test_bootstrap_scripts_use_normal_build_and_fast_suite(), test_clean_scripts_announce_repair_only_scope(), test_docs_make_confidence_suite_the_default_runtime_check(), test_github_default_events_do_not_run_duplicate_heavy_smokes(), test_github_default_smoke_uses_downstream_path_install_not_wheel_build(), test_github_full_packaging_remains_manual_only()

### Community 69 - "Community 69"
Cohesion: 0.33
Nodes (5): Return the dielectric model evaluation for the current state., One relative-permittivity data residual for regression/provenance workflows., Return a flat regression record with ``x_*`` composition columns., Return this residual as a first-class regression term descriptor., relative_permittivity()

### Community 70 - "Community 70"
Cohesion: 0.52
Nodes (6): _salt_mixture(), test_electrolyte_bubble_pressure_best_effort_returns_diagnostics(), test_electrolyte_bubble_pressure_converges_for_water_salt(), test_electrolyte_bubble_pressure_nonconvergence_raises_by_default(), test_electrolyte_bubble_pressure_rejects_ionic_vapor_species(), test_electrolyte_bubble_pressure_rejects_python_backend_alias()

### Community 71 - "Community 71"
Cohesion: 0.4
Nodes (6): il_mole_fraction(), scan_temperature_branch(), solve_binary_lle(), split_branches(), water_mole_fraction(), water_solubility_in_il()

### Community 72 - "Community 72"
Cohesion: 0.33
Nodes (5): Create an immutable thermodynamic state for the mixture.          States built, Return pressure-consistency diagnostics for an externally supplied density., Return the pressure of the bound state., pressure(), state()

### Community 73 - "Community 73"
Cohesion: 0.47
Nodes (3): test_generated_output_roots_are_not_tracked_in_analyses(), test_old_gallery_and_script_roots_are_not_tracked(), _tracked_files()

### Community 74 - "Community 74"
Cohesion: 0.7
Nodes (4): test_ternary_hydrocarbon_basis_tp_flash_closes_material_and_fugacity_balance(), test_tp_flash_phase_diagnostics_are_json_serializable_when_requested(), test_tp_flash_reports_no_split_when_rachford_rice_has_no_bracket(), hydrocarbon_basis_mixture()

### Community 76 - "Community 76"
Cohesion: 0.4
Nodes (5): _default_species_entry(), molality_to_molefraction(), molefraction_to_molality(), Convert salt molality (mol/kg solvent) to species mole-fraction vector., Convert mole fractions to molality for 1:1 salt systems.

### Community 77 - "Community 77"
Cohesion: 0.7
Nodes (4): _load_backend(), test_pep517_build_backend_honors_persistent_build_dir_env(), test_pep517_build_backend_preserves_explicit_build_dir(), test_pep517_build_backend_uses_isolated_build_dir_by_default()

### Community 78 - "Community 78"
Cohesion: 0.7
Nodes (4): _load_script(), test_build_script_auto_prefers_ninja_for_new_build_tree(), test_build_script_preserves_existing_generator_for_auto(), test_build_script_rejects_explicit_generator_switch_without_clean()

### Community 79 - "Community 79"
Cohesion: 0.83
Nodes (3): _draw_grid(), _to_xy(), _plot_tie_lines()

## Knowledge Gaps
- **191 isolated node(s):** `Scan alternative Figure 6b accounting methods against digitized contribution lin`, `Generate Figure 6b fit-check plots from bookkeeping curves plus z-correction wei`, `Diagnose exact Figure 6b bookkeeping, including the compressibility correction.`, `Generate Figure 6b fit-comparison plots on the paper-style ``mu`` basis.  Each`, `Scan Figure 6b contribution sensitivity to imposed liquid-density scaling.` (+186 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **15 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `ValueError()` connect `Community 0` to `Community 1`, `Community 2`, `Community 3`, `Community 4`, `Community 5`, `Community 6`, `Community 7`, `Community 8`, `Community 9`, `Community 10`, `Community 11`, `Community 12`, `Community 15`, `Community 16`, `Community 17`, `Community 21`, `Community 24`, `Community 28`, `Community 29`, `Community 30`, `Community 31`, `Community 33`, `Community 38`, `Community 39`, `Community 41`, `Community 42`, `Community 44`, `Community 46`, `Community 47`, `Community 54`, `Community 61`, `Community 72`, `Community 76`?**
  _High betweenness centrality (0.344) - this node is a cross-community bridge._
- **Why does `main()` connect `Community 9` to `Community 0`, `Community 2`, `Community 3`, `Community 6`, `Community 7`, `Community 8`, `Community 11`, `Community 12`, `Community 15`, `Community 17`, `Community 27`, `Community 29`, `Community 30`, `Community 31`, `Community 33`, `Community 35`, `Community 37`, `Community 38`, `Community 42`, `Community 43`, `Community 44`, `Community 55`, `Community 71`, `Community 79`?**
  _High betweenness centrality (0.247) - this node is a cross-community bridge._
- **Why does `InputError` connect `Community 28` to `Community 5`, `Community 7`, `Community 13`, `Community 16`, `Community 18`, `Community 19`, `Community 21`, `Community 23`, `Community 24`, `Community 26`, `Community 37`, `Community 41`, `Community 45`, `Community 47`, `Community 52`, `Community 53`, `Community 54`, `Community 57`, `Community 66`, `Community 69`, `Community 72`?**
  _High betweenness centrality (0.118) - this node is a cross-community bridge._
- **Are the 9 inferred relationships involving `main()` (e.g. with `il_label()` and `scan_temperature_branch()`) actually correct?**
  _`main()` has 9 INFERRED edges - model-reasoned connections that need verification._
- **Are the 158 inferred relationships involving `ValueError()` (e.g. with `species_for_combo()` and `normalized_comp()`) actually correct?**
  _`ValueError()` has 158 INFERRED edges - model-reasoned connections that need verification._
- **Are the 26 inferred relationships involving `InputError` (e.g. with `electrolyte_bubble.py` and `ElectrolyteBubbleResult`) actually correct?**
  _`InputError` has 26 INFERRED edges - model-reasoned connections that need verification._
- **Are the 5 inferred relationships involving `ePCSAFTState` (e.g. with `ActivityCoefficientResult` and `InputError`) actually correct?**
  _`ePCSAFTState` has 5 INFERRED edges - model-reasoned connections that need verification._