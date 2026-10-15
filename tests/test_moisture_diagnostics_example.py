"""Moisture, diagnostics and the worked example of the software paper."""
import os
import runpy
import math
import pytest
from rothcruz import dead_fuel_moisture, required_crown_fraction, simard_emc, wind_multiplier

EXAMPLE = os.path.join(os.path.dirname(__file__), "..", "examples", "uiseong_2025_interval.py")


def test_simard_branches_and_offsets():
    # the three regression branches of Simard (1968), evaluated by hand at 68 F (20 C)
    assert float(simard_emc(5.0, 20.0)) == pytest.approx(0.03229 + 0.281073 * 5 - 0.000578 * 5 * 68)
    assert float(simard_emc(30.0, 20.0)) == pytest.approx(2.22749 + 0.160107 * 30 - 0.01478 * 68)
    assert float(simard_emc(80.0, 20.0)) == pytest.approx(21.0606 + 0.005565 * 80 ** 2 - 0.00035 * 80 * 68 - 0.483199 * 80)
    m1, m10, m100 = dead_fuel_moisture(30.0, 20.0)
    assert float(m1) == pytest.approx(1.03 * float(simard_emc(30.0, 20.0))) and float(m10 - m1) == pytest.approx(1) and float(m100 - m1) == pytest.approx(2)


def test_diagnostics_basic():
    assert required_crown_fraction(10.0, [1.0, 1.0], [11.0, 11.0]) == pytest.approx(0.9)
    assert wind_multiplier(8.0, lambda z: 2.0 * z) == pytest.approx(4.0, rel=1e-8)
    assert math.isnan(wind_multiplier(1.0, lambda z: 2.0 + z))         # model already above the observation at lo


def test_worked_example_reproduces_the_paper():
    out = runpy.run_path(EXAMPLE)["main"](verbose=False)
    assert out["n"] == 21
    assert out["surface"] == pytest.approx(1.164, abs=5e-3) and out["r91"] == pytest.approx(10.145, abs=5e-3)
    assert out["c05"] == pytest.approx(50.03, abs=5e-2) and out["scott_reinhardt"] == pytest.approx(6.05, abs=5e-3)
    assert out["f_r91"] == pytest.approx(7.23, abs=5e-3) and out["f_c05"] == pytest.approx(1.33, abs=5e-3)
    assert out["m_s"] == pytest.approx(51.4, abs=0.05)
    assert (out["torching"], out["crowning"]) == (12, 2) and out["cfb"] == pytest.approx(0.32, abs=5e-3)
