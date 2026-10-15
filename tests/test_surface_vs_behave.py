"""Surface fire against the BehavePlus engine (firelab/behave, commit d4acc8a); see tests/reference/README.md."""
import os
import numpy as np
import pandas as pd
from rothcruz import BEHAVE_COMPAT, PUBLISHED, fuel_model, rothermel

REF = os.path.join(os.path.dirname(__file__), "reference")


def _rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = np.abs(a - b) / np.maximum(np.abs(b), 1e-12)
    d[(np.abs(a) < 1e-12) & (np.abs(b) < 1e-12)] = 0
    return d.max()


def test_surface_vs_behaveplus():
    ref = pd.read_csv(os.path.join(REF, "behave_surface_reference.csv"))
    assert len(ref) == 7632 and ref.fm.nunique() == 53
    worst = 0.0
    for fm, g in ref.groupby("fm"):
        r = rothermel(fuel_model(fm, BEHAVE_COMPAT), g.m1 / 100, g.m10 / 100, g.m100 / 100,
                      g.mlh / 100, g.mlw / 100, g.wind_mph * 88.0, g.slope_pct / 100, options=BEHAVE_COMPAT)
        for k, col in [("ros", "ros_ftmin"), ("ib", "ib_btufts"), ("ir", "ir_btuft2min"), ("flame", "fl_ft")]:
            worst = max(worst, _rel(r[k], g[col]))
    assert worst < 1e-5, worst


def test_steep_slope_cap_reproduces_engine():
    """D7: the engine caps the slope factor at 0.9 * reaction intensity; BEHAVE_COMPAT reproduces it, PUBLISHED does not."""
    ref = pd.read_csv(os.path.join(REF, "behave_steep_slope_reference.csv")).rename(columns={"limit_0.9IR": "limit_0_9IR"})
    capped = 0
    for fm, g in ref.groupby("fuel_model"):
        f = fuel_model(int(fm), BEHAVE_COMPAT)
        args = (.06, .07, .08, .6, .9)
        r0 = rothermel(f, *args, 0.0, 0.0, options=BEHAVE_COMPAT)["ros"]
        for row in g.itertuples():
            tan_s = row.slope_pct / 100.0
            phi_c = float(rothermel(f, *args, 0.0, tan_s, options=BEHAVE_COMPAT)["ros"] / r0 - 1)
            phi_p = float(rothermel(f, *args, 0.0, tan_s, options=PUBLISHED)["ros"] / r0 - 1)
            assert phi_c == np.float64(row.phi_s_engine) or abs(phi_c / row.phi_s_engine - 1) < 1e-3, (fm, row.slope_pct)
            if row.phi_s_predicted_by_tan2 > row.limit_0_9IR * 1.001:      # the cap is active in the engine
                capped += 1
                assert phi_p > phi_c * 1.001                               # published value is not capped
    assert capped >= 6


def test_exact_curing_differs_only_for_dynamic_models():
    """Documents the BehavePlus rounded curing formula (1.333-1.11*M), divergence D1."""
    rounded = PUBLISHED.replace(curing="behaveplus")
    for fm in (1, 10, 165):          # static models: identical
        a = rothermel(fuel_model(fm), .06, .07, .08, .6, .9, 440, 0, options=PUBLISHED)["ros"]
        b = rothermel(fuel_model(fm), .06, .07, .08, .6, .9, 440, 0, options=rounded)["ros"]
        assert abs(a - b) < 1e-12
    # GR6 at Mlh = 120 %: exact formula gives 0 cured, BehavePlus 0.001 cured.
    # The tiny transferred load matters because GR6 has almost no dead fuel.
    a = rothermel(fuel_model(106), .06, .07, .08, 1.2, 1.5, 440, 0, options=PUBLISHED)["ib"]
    b = rothermel(fuel_model(106), .06, .07, .08, 1.2, 1.5, 440, 0, options=rounded)["ib"]
    assert abs(a - b) / b > 1e-2
