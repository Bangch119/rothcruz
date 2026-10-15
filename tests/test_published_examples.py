"""Worked examples of Cruz et al. (2005)."""
import numpy as np
import pytest


def test_cruz2005_published_examples():
    from rothcruz import cruz2005_crown_ros
    # Cruz et al. (2005) text: EFFM 4 %, CBD 0.2, U10 50 and 75 km/h -> 137 and 196 m/min
    # The published coefficients are rounded (11.02, 0.90, 0.19, 0.17); with them the paper's
    # own worked examples are reproduced to within 3 % (all ~1.5-2.5 % high).
    r = cruz2005_crown_ros([50.0, 75.0], 0.2, 4.0)
    assert np.allclose(r["active"], [137.0, 196.0], rtol=0.03)
    assert float(cruz2005_crown_ros(7.0, 0.08, 8.0)["passive"]) == pytest.approx(7.6, rel=0.03)    # Bor Island
    assert float(cruz2005_crown_ros(21.0, 0.1, 10.0)["passive"]) == pytest.approx(10.2, rel=0.03)  # Ocala
    assert np.allclose(cruz2005_crown_ros(74.0, [0.05, 0.1], 9.0)["active"], [63.2, 72.2], rtol=0.03)  # Burnt fire

