from __future__ import annotations

import epcsaft


def test_runtime_reports_cppad_build_contract() -> None:
    info = epcsaft.runtime_build_info()
    cppad = info["optional_dependencies"]["cppad"]

    assert cppad["backend"] == "cppad"
    assert cppad["status"] == "enabled_available"
    assert cppad["compiled"] is True
    assert cppad["available"] is True

    capabilities = epcsaft.capabilities()
    assert capabilities["derivatives"]["cppad"]["status"] == cppad["status"]
    assert capabilities["derivatives"]["cppad"]["compiled"] is cppad["compiled"]
    assert capabilities["derivatives"]["cppad"]["production"] is False
    assert capabilities["derivatives"]["cppad"]["production_eos_coverage"] is False
