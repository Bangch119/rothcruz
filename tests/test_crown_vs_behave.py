"""Crown fire against the BehavePlus engine (firelab/behave, commit d4acc8a)."""
import os
import numpy as np
import pandas as pd
import pytest
from rothcruz import BEHAVE_COMPAT, PUBLISHED, crown_fire, critical_active_ros, critical_intensity

REF = os.path.join(os.path.dirname(__file__), "reference")


def _rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = np.abs(a - b) / np.maximum(np.abs(b), 1e-12)
    d[(np.abs(a) < 1e-12) & (np.abs(b) < 1e-12)] = 0
    return d.max()


def test_crown_vs_behaveplus():
    ref = pd.read_csv(os.path.join(REF, "behave_crown_reference.csv"))
    mism, w_s, w_c, w_i = 0, 0, 0, 0
    for fm, g in ref.groupby("fm"):
        r = crown_fire(fm, g.m1, g.m10, g.m100, g.mlh, g.mlw, g.fmc, g.wind20_kmh, g.slope_pct,
                       g.cbh_m, g.cbd, g.cc_pct, g.ch_m, None, options=BEHAVE_COMPAT)
        w_s = max(w_s, _rel(r["surface_ros"], g.surf_ros_mmin))
        w_c = max(w_c, _rel(r["crown_ros"], g.crown_ros_mmin))
        w_i = max(w_i, _rel(r["crown_I"], g.crown_ib_kwm))
        mism += int((r["fire_type"] != g.fire_type.values).sum())
    assert set(ref.fire_type) == {0, 1, 2, 3}       # all four regimes exercised
    assert w_s < 1e-5 and w_c < 1e-6 and w_i < 1e-4 and mism == 0


def test_van_wagner_analytic():
    # CBH 4 m, FMC 100 % : (0.010*4*(460+2590))**1.5 = 122**1.5
    assert critical_intensity(4.0, 100.0) == pytest.approx(122.0 ** 1.5, rel=1e-12)
    assert critical_active_ros(0.15) == pytest.approx(20.0)
    # D5: the engine limits CBH >= 0.1 m and FMC >= 30 %; the published equation has no limits
    assert critical_intensity(0.0, 10.0, limits=True) == pytest.approx((0.001 * (460 + 25.9 * 30)) ** 1.5)
    assert critical_intensity(0.0, 10.0) == 0.0
