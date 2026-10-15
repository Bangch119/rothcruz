"""The options object: every divergence from the engine can be switched, and the presets differ where documented."""
import numpy as np
import pytest
from rothcruz import (BEHAVE_COMPAT, FUEL_MODELS, PUBLISHED, ModelOptions, crown_fire, fuel_model, rothermel,
                      scott_reinhardt, surface_fire)
from rothcruz.options import resolve

ARGS = (.06, .07, .08, .6, .9)         # dead moistures, live herbaceous, live woody


def test_presets_and_validation():
    assert resolve(None) is PUBLISHED and resolve("published") is PUBLISHED and resolve("behave_compat") is BEHAVE_COMPAT
    assert PUBLISHED.curing == "exact" and PUBLISHED.wind_limit == "none" and PUBLISHED.fuel_table == "published"
    assert BEHAVE_COMPAT.curing == "behaveplus" and BEHAVE_COMPAT.fuel_table == "engine" and BEHAVE_COMPAT.engine_constants
    assert PUBLISHED.replace(wind_limit="effective").wind_limit == "effective"
    for bad in (dict(curing="x"), dict(wind_limit="andrews"), dict(fuel_table="y")):
        with pytest.raises(ValueError):
            ModelOptions(**bad)
    with pytest.raises(ValueError):
        resolve("nonsense")


def test_published_fuel_table_differs_from_the_engine_only_for_fm7():
    """D9: Albini (1976) Table 7 gives a live SAV of 1550 for model 7; the engine has 1500."""
    differ = [n for n in range(1, 205) if n in FUEL_MODELS and fuel_model(n, PUBLISHED) != fuel_model(n, BEHAVE_COMPAT)]
    assert differ == [7]
    assert fuel_model(7, PUBLISHED).savwoody == 1550 and fuel_model(7, BEHAVE_COMPAT).savwoody == 1500
    a = float(rothermel(fuel_model(7, PUBLISHED), *ARGS, 440.0, 0.0)["ros"])
    b = float(rothermel(fuel_model(7, BEHAVE_COMPAT), *ARGS, 440.0, 0.0)["ros"])
    assert 0.990 < a / b < 0.999            # 1550 lowers the spread rate by about 0.7 % (0.5-1.5 % over the grid)
    custom = fuel_model(7, BEHAVE_COMPAT)   # a model the user passes in is never replaced
    assert fuel_model(custom, PUBLISHED) is custom


@pytest.mark.parametrize("model", [1, 10, 165])
def test_wind_limit_options(model):
    fm = fuel_model(model)
    ir = float(rothermel(fm, *ARGS, 440.0, 0.0)["ir"])
    u_lim = 0.9 * ir
    eff, mid = PUBLISHED.replace(wind_limit="effective"), PUBLISHED.replace(wind_limit="midflame")
    # below the limit nothing changes
    base = float(rothermel(fm, *ARGS, 0.5 * u_lim, 0.0)["ros"])
    assert float(rothermel(fm, *ARGS, 0.5 * u_lim, 0.0, options=eff)["ros"]) == pytest.approx(base, rel=1e-12)
    # far above the limit the spread rate saturates at the value reached at the limit
    at_limit = float(rothermel(fm, *ARGS, u_lim, 0.0, options=eff)["ros"])
    far = rothermel(fm, *ARGS, 100.0 * u_lim, 0.0, options=eff)
    assert float(far["ros"]) == pytest.approx(at_limit, rel=1e-12) and bool(far["wind_limited"])
    assert float(rothermel(fm, *ARGS, 100.0 * u_lim, 0.0)["ros"]) > 5 * at_limit      # no limit: keeps growing
    # on flat ground limiting the effective wind and the midflame wind are the same
    assert float(rothermel(fm, *ARGS, 100.0 * u_lim, 0.0, options=mid)["ros"]) == pytest.approx(at_limit, rel=1e-12)
    # with slope, the effective-wind limit is the stricter one, and spread rate never decreases with wind
    s_eff = [float(rothermel(fm, *ARGS, u, 0.5, options=eff)["ros"]) for u in np.linspace(0, 20 * u_lim, 40)]
    assert np.all(np.diff(s_eff) >= -1e-9)
    assert (float(rothermel(fm, *ARGS, 5 * u_lim, 0.5, options=eff)["ros"])
            < float(rothermel(fm, *ARGS, 5 * u_lim, 0.5, options=mid)["ros"])
            <= float(rothermel(fm, *ARGS, 5 * u_lim, 0.5)["ros"]))


def test_crown_ratio_override_and_default_wind_limit_off():
    kw = dict(fuel=165, m1=4, m10=5, m100=6, mlh=30, mlw=90, fmc=100, wind20_kmh=30.0, slope_pct=20.0, cbh_m=2.0,
              cbd_kgm3=0.15, canopy_cover_pct=50, canopy_height_m=15.0)
    pub_user = crown_fire(**kw, crown_ratio=0.2, options=PUBLISHED)["waf"]
    pub_derived = crown_fire(**kw, crown_ratio=None, options=PUBLISHED)["waf"]
    com_user = crown_fire(**kw, crown_ratio=0.2, options=BEHAVE_COMPAT)["waf"]
    com_derived = crown_fire(**kw, crown_ratio=None, options=BEHAVE_COMPAT)["waf"]
    assert pub_user != pub_derived                      # published: the user's crown ratio is used
    assert com_user == com_derived == pub_derived       # engine: it is overridden by (CH - CBH)/CH


def test_scott_reinhardt_modes_are_consistent():
    kw = dict(fuel=165, m1=4, m10=5, m100=6, mlh=30, mlw=90, fmc=100, wind20_kmh=np.array([20.0, 30.0, 44.0]),
              slope_pct=20.0, cbh_m=2.0, cbd_kgm3=0.15, canopy_cover_pct=50, canopy_height_m=15.0)
    p, c = scott_reinhardt(**kw, options=PUBLISHED), scott_reinhardt(**kw, options=BEHAVE_COMPAT)
    assert np.all((p["cfb"] >= 0) & (p["cfb"] <= 1)) and np.all((c["cfb"] >= 0) & (c["cfb"] <= 1))
    assert np.allclose(p["final_ros"], c["final_ros"], rtol=1e-3)        # the two modes agree closely for this fuel model
    out = surface_fire(165, 4, 5, 6, 30, 90, 2.0, 10.0, options=BEHAVE_COMPAT)
    assert out["ros"] > 0
