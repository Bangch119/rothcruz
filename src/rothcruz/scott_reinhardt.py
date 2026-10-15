"""
Scott & Reinhardt (2001) linkage of surface and crown fire models: crown fraction burned (CFB),
critical open wind speed O'active, final spread rate and fireline intensity.

Calculation order follows the BehavePlus engine (firelab/behave, crown.cpp, doCrownRunScottAndReinhardt):
  1. surface run (site WAF, slope)                      -> Rs, HPA_s, I_s
  2. fuel model 10 run (WAF 0.4, no slope)               -> R_active = 3.34 * R10
  3. Van Wagner I0 -> R'initiation = 60 * I0 / HPA_s
  4. R'active = 3 / CBD -> O'active (20-ft wind at which 3.34*R10 = R'active)
  5. R'SA = surface spread rate at wind O'active
  6. CFB = (Rs - R'initiation) / (R'SA - R'initiation), limited to [0, 1]
  7. final: surface/conditional -> Rs; torching -> Rs + CFB*(R_active - Rs); crowning -> R_active
options.engine_constants (BEHAVE_COMPAT) reproduces the engine's rounded constants (16.0185, 3.28084, 0.2886719, 3.464140419)
and its classification order; the other options are those of crown_fire() and rothermel().
"""
from __future__ import annotations
import numpy as np
from .constants import (BTU_FT_S_TO_KW_M_ENGINE, FT_PER_M_ENGINE, KGM3_TO_LBFT3_ENGINE, KW_M_TO_BTU_FT_S_ENGINE)
from .options import resolve
from .surface import BTU_FT_S_TO_KW_M, FT_PER_M, fuel_model, rothermel
from .crown import KGM3_TO_LBFT3, KMH_TO_FTMIN, critical_intensity, crown_fire


def scott_reinhardt(fuel, m1, m10, m100, mlh, mlw, fmc, wind20_kmh, slope_pct, cbh_m, cbd_kgm3,
                    canopy_cover_pct, canopy_height_m, crown_ratio=None, *, options=None):
    """
    All inputs scalars or arrays (moistures in %, 20-ft open wind in km/h). Returns a dict, SI units.
    options: ModelOptions, "published", "behave_compat" or None (= published). The wind-limit option applies to every
    Rothermel surface calculation in the chain, including the fuel model 10 runs.
    """
    o = resolve(options)
    ec = o.engine_constants
    fm = fuel_model(fuel, o); fm10 = fuel_model(10, o)
    pct = lambda x: np.asarray(x, float) / 100.0
    base = crown_fire(fuel, m1, m10, m100, mlh, mlw, fmc, wind20_kmh, slope_pct, cbh_m, cbd_kgm3,
                      canopy_cover_pct, canopy_height_m, crown_ratio, options=o)
    U20 = np.asarray(wind20_kmh, float) * KMH_TO_FTMIN
    waf = base["waf"]; M = [pct(x) for x in (m1, m10, m100, mlh, mlw)]; tanS = pct(slope_pct)
    surf = rothermel(fm, *M, waf * U20, tanS, options=o)
    rs, hpa_s = surf["ros"], surf["hpa"]                                  # ft/min, Btu/ft2
    tan10 = tanS if o.crown_slope else 0.0                                     # D10
    r_active = 3.34 * rothermel(fm10, *M, 0.4 * U20, tan10, options=o)["ros"]   # ft/min

    # Van Wagner critical intensity and Scott & Reinhardt R'initiation
    i0 = critical_intensity(cbh_m, fmc, limits=o.van_wagner_limits)         # kW/m
    i0_btu = i0 * (KW_M_TO_BTU_FT_S_ENGINE if ec else 1.0 / BTU_FT_S_TO_KW_M)
    r_init = np.where(hpa_s > 0, 60.0 * i0_btu / np.maximum(hpa_s, 1e-12), 0.0)

    # R'active (ft/min) and O'active from the fuel model 10 wind function
    cbd = np.asarray(cbd_kgm3, float)
    kg2lb = KGM3_TO_LBFT3_ENGINE if ec else KGM3_TO_LBFT3
    cbd_eff = (16.0185 if ec else 16.01846337) * cbd * kg2lb                # kg/m3 as the engine sees it
    r_act_crit = (FT_PER_M_ENGINE if ec else FT_PER_M) * 3.0 / np.maximum(cbd_eff, 1e-12)
    r10 = r_act_crit / 3.34
    f10 = rothermel(fm10, *M, 0.0, 0.0, options=o)
    ros0 = f10["ir"] * f10["xi"] / f10["sink"]
    K = f10["wind_C"] * np.power(f10["rel_packing"], -f10["wind_E"])
    phi_s10 = rothermel(fm10, *M, 0.0, tan10, options=o)["phi_s"] if o.crown_slope else 0.0
    a = (r10 / ros0 - 1.0 - phi_s10) / K                                       # Scott & Reinhardt (2001), Eq. 19
    o_mid = np.where(a > 0, np.power(np.maximum(a, 0.0), 1.0 / f10["wind_B"]), 0.0)   # slope alone sufficient: O'active = 0
    o_active = o_mid / 0.4                                                  # ft/min, 20-ft wind
    r_sa = rothermel(fm, *M, waf * o_active, tanS, options=o)["ros"]        # ft/min

    den = r_sa - r_init
    cfb = np.where(den > 1e-7, (rs - r_init) / np.where(den > 1e-7, den, 1.0), 0.0)
    cfb = np.clip(cfb, 0.0, 1.0)

    ft = np.asarray(base["fire_type"])          # classification of crown_fire(), same options
    ch_ft = np.asarray(canopy_height_m, float) * FT_PER_M
    canopy_hpa = cbd * kg2lb * (ch_ft - np.asarray(cbh_m, float) * FT_PER_M) * 8000.0
    passive_ros = rs + cfb * (r_active - rs)
    passive_hpa = hpa_s + canopy_hpa * cfb
    crown_hpa = hpa_s + canopy_hpa
    final_ros = np.where(ft == 3, r_active, np.where(ft == 1, passive_ros, rs))
    final_ib_btu = np.where(ft == 3, r_active / 60.0 * crown_hpa, np.where(ft == 1, passive_hpa * passive_ros / 60.0, surf["ib"]))
    kw = BTU_FT_S_TO_KW_M_ENGINE if ec else BTU_FT_S_TO_KW_M
    return dict(cfb=cfb, final_ros=final_ros / FT_PER_M, final_ib=final_ib_btu * kw, passive_ros=passive_ros / FT_PER_M,
                o_active_kmh=o_active / KMH_TO_FTMIN, r_initiation=r_init / FT_PER_M, r_sa=r_sa / FT_PER_M,
                surface_ros=rs / FT_PER_M, crown_ros=r_active / FT_PER_M, fire_type=ft)
