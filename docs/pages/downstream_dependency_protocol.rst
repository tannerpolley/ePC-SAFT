Downstream Dependency Protocol
==============================

Use this protocol when another project fails because of a suspected upstream
``epcsaft`` package bug, missing public API contract, packaging problem,
documentation ambiguity, or numerical behavior change.

Issue intake
------------

Downstream projects should file GitHub issues against ``tannerpolley/ePC-SAFT``
with the downstream dependency bug or upstream package request form. A complete
report includes the downstream repo or path, task goal, failing command, error
or bad result, minimal reproducer, imported ``epcsaft`` path, expected behavior,
actual behavior, downstream validation command after the fix, and temporary
workaround status.

Use labels to make routing durable:

``downstream-bug``
   The issue came from a downstream project and may require an upstream package
   fix.

``needs-repro``
   The report is missing enough evidence to reproduce or validate the package
   behavior.

``agent-ready``
   The report has a minimal reproducer, imported package path, and downstream
   validation command.

``in-progress``, ``upstream-fix-ready``, ``downstream-validated``, and
``blocked-downstream`` track the handoff state. Area labels such as
``python-api``, ``native``, ``solver``, ``packaging``, ``docs``, ``validation``,
``regression``, and ``equilibrium`` identify the likely owner.

Upstream workflow
-----------------

For complete reports, reproduce the failure before editing. Classify the issue
as API misuse, missing public contract, native equation issue, packaging/build
issue, numerical/solver issue, or documentation ambiguity. Add a focused
regression test when feasible, fix the package rather than patching downstream
code, and validate with the smallest relevant ladder:

.. code-block:: powershell

   uv run python run_pytest.py <focused-test-targets> -q
   uv run python scripts/validate_project.py quick

Run ``uv run python scripts/build_epcsaft.py`` and ``uv run python
scripts/doctor.py`` when native code or package import state is involved. Use
``uv run python scripts/validate_project.py confidence`` only for broad or
native-risk changes.

Before handing the issue back, comment with the root cause, files or APIs
changed, focused tests run, install or rebuild command, downstream validation
command, and any remaining risk. Mark the issue ``upstream-fix-ready`` until
the downstream project confirms the original command now passes, then mark it
``downstream-validated`` and close it.

Local triage helper
-------------------

Use the read-only helper to inspect one issue and choose the first command:

.. code-block:: powershell

   uv run python scripts/triage_dependency_issue.py --issue 12
   uv run python scripts/triage_dependency_issue.py --issue https://github.com/tannerpolley/ePC-SAFT/issues/12 --json

The helper fetches the issue through GitHub CLI, checks the required sections,
classifies the likely area from labels and body text, and prints recommended
next commands. It exits nonzero when the report still needs a reproducer,
imported package path, or downstream validation command.

Automation
----------

A local Codex heartbeat can poll open ``downstream-bug`` and ``agent-ready``
issues every five minutes. It should ignore issues labeled ``in-progress``,
``blocked-downstream``, or ``downstream-validated``. For ``needs-repro`` issues,
it should report missing fields and avoid code changes. For ``agent-ready``
issues, it should run the triage helper first, then start the normal upstream
fix workflow only when the report is complete.
