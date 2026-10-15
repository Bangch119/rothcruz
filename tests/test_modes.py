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
    p, c = scott_reinhardt(**kw, options=PUBLISHED.replace(crown_slope=False)), scott_reinhardt(**kw, options=BEHAVE_COMPAT)
    assert np.all((p["cfb"] >= 0) & (p["cfb"] <= 1)) and np.all((c["cfb"] >= 0) & (c["cfb"] <= 1))
    assert np.allclose(p["final_ros"], c["final_ros"], rtol=1e-3)        # the two modes agree closely for this fuel model
    out = surface_fire(165, 4, 5, 6, 30, 90, 2.0, 10.0, options=BEHAVE_COMPAT)
    assert out["ros"] > 0


def test_crown_slope_d10():
    """D10: Rothermel (1991) and Scott & Reinhardt (2001, Eq. 8, 19) include the slope factor of fuel model 10 in the crown
    fire spread rate; the engine (BEHAVE_COMPAT) uses zero slope."""
    from rothcruz import rothermel, fuel_model
    from rothcruz.crown import KMH_TO_FTMIN
    kw = dict(fuel=165, m1=4, m10=5, m100=6, mlh=30, mlw=90, fmc=100, wind20_kmh=25.0, cbh_m=2.0, cbd_kgm3=0.15,
              canopy_cover_pct=50, canopy_height_m=15.0)
    flat = crown_fire(**kw, slope_pct=0.0)["crown_ros"]
    steep = crown_fire(**kw, slope_pct=40.0)["crown_ros"]
    r = rothermel(fuel_model(10), .04, .05, .06, .30, .90, 0.4 * 25.0 * KMH_TO_FTMIN, 0.40)
    assert float(steep / flat) == pytest.approx(float((1 + r["phi_w"] + r["phi_s"]) / (1 + r["phi_w"])), rel=1e-9)
    for opts in (BEHAVE_COMPAT, PUBLISHED.replace(crown_slope=False)):          # slope ignored in the crown rate
        assert float(crown_fire(**kw, slope_pct=40.0, options=opts)["crown_ros"]) == pytest.approx(float(flat), rel=1e-12)
    # O'active is consistent: at the 20-ft wind O'active the crown rate (with slope) equals R'active = 3/CBD
    for slope in (0.0, 30.0):
        o = scott_reinhardt(**kw, slope_pct=slope)["o_active_kmh"]
        at = crown_fire(**{**kw, "wind20_kmh": float(o)}, slope_pct=slope)
        assert float(at["crown_ros"]) == pytest.approx(3.0 / 0.15, rel=1e-6)
    assert float(scott_reinhardt(**kw, slope_pct=30.0)["o_active_kmh"]) < float(scott_reinhardt(**kw, slope_pct=0.0)["o_active_kmh"])
