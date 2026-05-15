# Native EOS, Association Sensitivities, And Ceres Benchmark Goal

Source authority: GitHub issue #114, "Complete native EOS/AD, association implicit sensitivities, and associating-binary Ceres regression"

Issue URL: https://github.com/tannerpolley/ePC-SAFT/issues/114

## Outcome

Execute issue #114 as a vertical implementation goal. The package must support one complete production chain: native EOS contribution derivatives, association site-fraction implicit sensitivities, associating-mixture binary-interaction total derivatives, native Ceres regression, and an associating-binary benchmark proof through the public Python regression API.

This board setup does not complete issue #114. It creates the execution board and the issue-required intake note path so a later `/goal` run can proceed through the issue phases in order.

## Hard Constraints

- Follow issue #114 exactly. If this file and the issue disagree, the GitHub issue is authoritative.
- Do not start source edits until `docs/goals/native-eos-association-ceres-benchmark/notes/intake.md` is fully populated with the Stage 0 evidence required by issue #114.
- Do not close, mark ready, or call issue #114 complete from inventories, manifests, schema-only work, diagnostics-only routes, staged-only helpers, synthetic payloads, mocked solver output, documented limitations, or unsupported capability labels.
- Use only production derivative mechanisms allowed by issue #114: analytic, CppAD, analytic implicit sensitivity, or CppAD implicit sensitivity.
- Ceres must own the production optimizer loop for the required regression proof. Python may validate inputs, call native code, and format outputs, but must not own the production objective loop.
- Keep `epcsaft` a general-purpose package. Do not add public APIs named after downstream applications or downstream metrics.
- Do not edit equilibrium solver implementation except for read-only API compatibility fixes explicitly justified by the issue work.
- Avoid committing the old missing-backend token or old numeric-derivative wording as contiguous literals. Construct guard searches from fragments if needed.

## Required Phase Order

1. Stage 0 - Intake gate: fill `notes/intake.md` with the branch/base snapshot, inspected native EOS files, inspected association solver files, inspected native regression files, current failing or skipped associating-binary regression path, and chosen benchmark fixture.
2. Pre-Stage 1 hard requirement - remove the legacy Eigen AD production route and move the required derivative path to CppAD or implicit sensitivities.
3. Stage 1 - Native EOS/AD substrate: harden scalar-templated native evaluation and derivative result shape for the benchmark path.
4. Stage 2 - Association site-fraction implicit sensitivities: implement total derivatives through solved association site fractions.
5. Stage 3 - Binary interaction derivative propagation: prove associating-mixture `k_ij` derivative participation and stable Jacobian shape; handle `l_ij` and `k_hb_ij` according to issue #114.
6. Stage 4 - Native Ceres regression: make native Ceres own the production fit loop and return the required fit result fields.
7. Stage 5 - Benchmark proof: run an associating-binary `k_ij` regression benchmark using repo-contained data and write `notes/final_benchmark_report.md`.
8. Required validation - run the exact build, pytest, docs, quick-validation, and whitespace checks named in issue #114.
9. Final audit - verify every issue #114 completion line is true before calling the issue complete; otherwise leave the issue open or PR draft with a stopped-state report.

## Completion Proof

The goal is complete only when every issue #114 stage is implemented, the associating-binary benchmark proves the native Ceres plus CppAD/implicit derivative chain, required tests and docs pass, public Python regression API works, docs state exactly what is proven, and the final audit confirms no staged-only, diagnostic-only, manifest-only, synthetic-only, or limitation-only result was used.

## Starter Command

`/goal Follow docs/goals/native-eos-association-ceres-benchmark/goal.md.`
