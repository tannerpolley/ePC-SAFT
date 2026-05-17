from __future__ import annotations

import epcsaft


def test_runtime_reports_cppad_build_contract() -> None:
    info = epcsaft.runtime_build_info()
    cppad = info["native_dependencies"]["cppad"]

    assert cppad["backend"] == "cppad"
    assert cppad["status"] == "enabled_available"
    assert cppad["required"] is True
    assert cppad["compiled"] is True
    assert cppad["available"] is True

    capabilities = epcsaft.capabilities()
    assert capabilities["derivatives"]["cppad"]["status"] == cppad["status"]
    assert capabilities["derivatives"]["cppad"]["compiled"] is cppad["compiled"]
    assert "production" not in capabilities["derivatives"]["cppad"]
    assert "production_eos_coverage" not in capabilities["derivatives"]["cppad"]
