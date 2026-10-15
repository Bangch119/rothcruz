"""
Surface-to-crown transition and crown fire spread.

- Van Wagner (1977) critical surface fireline intensity for crown ignition
      I0 = (0.010 * CBH * (460 + 25.9 * FMC))**1.5     [kW/m]  (CBH m, FMC %)
- Van Wagner (1977) critical crown spread rate for active crowning
      R0 = 3.0 / CBD                                 [m/min] (CBD kg/m3)
- Rothermel (1991) crown spread rate
      R_crown = 3.34 * R_FM10(wind = 0.4 * U_20ft, slope = 0)
- Fire-type classification as in BehavePlus' Rothermel crown method
  (Surface / Torching / Conditional crown / Crowning)
- Wind adjustment factor, Albini & Baughman (1979) as given in Andrews (2012)
"""
from __future__ import annotations

import numpy as np

from .constants import (BTU_FT_S_TO_KW_M_ENGINE, FT_PER_M_ENGINE, KGM3_TO_LBFT3_ENGINE, KW_M_TO_BTU_FT_S_ENGINE)
from .options import resolve
from .surface import (FT_PER_M, BTU_FT_S_TO_KW_M, fuel_model, rothermel)

KMH_TO_FTMIN = 1000.0 * FT_PER_M / 60.0
KGM3_TO_LBFT3 = 0.0624279606
FIRE_TYPES = {0: "Surface", 1: "Torching", 2: "Conditional crown", 3: "Crowning"}


def critical_intensity(cbh_m, fmc_pct, *, limits=False):
    """Van Wagner (1977) I0 in kW/m. limits=True applies the engine's lower bounds (CBH 0.1 m, FMC 30 %), D5."""
    cbh = np.asarray(cbh_m, float)
    fmc = np.asarray(fmc_pct, float)
    if limits:
        cbh, fmc = np.maximum(cbh, 0.1), np.maximum(fmc, 30.0)
    return np.power(0.010 * cbh * (460.0 + 25.9 * fmc), 1.5)


def critical_active_ros(cbd_kgm3):
    """Van Wagner (1977) R0 = 3/CBD in m/min."""
    cbd = np.asarray(cbd_kgm3, float)
    return np.where(cbd > 1e-7, 3.0 / np.maximum(cbd, 1e-7), 0.0)


def wind_adjustment_factor(fuel_depth_ft, canopy_cover_frac, canopy_height_ft, crown_ratio):
    """Albini & Baughman (1979) WAF; sheltered if crown fill >= 0.05 and CH >= 6 ft."""
    cc = np.asarray(canopy_cover_frac, float)
    ch = np.asarray(canopy_height_ft, float)
    f = np.asarray(crown_ratio, float) * cc / 3.0
    sheltered = (cc >= 1e-7) & (f >= 0.05) & (ch >= 6.0)
    h = fuel_depth_ft
    unsh = 1.83 / np.log((20.0 + 0.36 * h) / (0.13 * h))
    chs = np.maximum(ch, 6.0)
    sh = 0.555 / (np.sqrt(np.maximum(f, 1e-9) * chs) * np.log((20.0 + 0.36 * chs) / (0.13 * chs)))
    return np.where(sheltered, sh, unsh)


def crown_fire(fuel, m1, m10, m100, mlh, mlw, fmc, wind20_kmh, slope_pct,
               cbh_m, cbd_kgm3, canopy_cover_pct, canopy_height_m, crown_ratio=None, *, options=None):
    """
    Combined surface + crown fire assessment (wind assumed upslope).
    Moistures in %, 20-ft (6.1 m) open wind in km/h.
    crown_ratio : a user value is used in published mode; if None it is derived as (CH - CBH)/CH. With
                  options.override_crown_ratio (engine) a user value is ignored, as BehavePlus does.
    options     : ModelOptions, "published", "behave_compat" or None (= published).
    Returns dict (SI): surface ros/I_B, I0, R0, crown ros, crown I_B,
    transition ratio, active ratio, fire_type code.
    """
    o = resolve(options)
    fm = fuel_model(fuel, o)
    pct = lambda x: np.asarray(x, float) / 100.0
    U20 = np.asarray(wind20_kmh, float) * KMH_TO_FTMIN
    ch_ft = np.asarray(canopy_height_m, float) * FT_PER_M
    if crown_ratio is None or o.override_crown_ratio:
        chm = np.asarray(canopy_height_m, float)
        crown_ratio = np.where(chm > 0, (chm - np.asarray(cbh_m, float)) / np.maximum(chm, 1e-9), 0.0)
    waf = wind_adjustment_factor(fm.depth, pct(canopy_cover_pct), ch_ft, crown_ratio)

    surf = rothermel(fm, pct(m1), pct(m10), pct(m100), pct(mlh), pct(mlw),
                     waf * U20, pct(slope_pct), options=o)
    fm10 = rothermel(fuel_model(10, o), pct(m1), pct(m10), pct(m100), pct(mlh), pct(mlw),
                     0.4 * U20, 0.0, options=o)
    crown_ros_ftmin = 3.34 * fm10["ros"]

    kw = BTU_FT_S_TO_KW_M_ENGINE if o.engine_constants else BTU_FT_S_TO_KW_M
    kg2lb = KGM3_TO_LBFT3_ENGINE if o.engine_constants else KGM3_TO_LBFT3
    I_s = surf["ib"] * kw
    I0 = critical_intensity(cbh_m, fmc, limits=o.van_wagner_limits)
    R0 = critical_active_ros(cbd_kgm3)
    crown_ros = crown_ros_ftmin / FT_PER_M

    canopy_load = (np.asarray(cbd_kgm3, float) * kg2lb
                   * (ch_ft - np.asarray(cbh_m, float) * FT_PER_M))
    crown_hpa = surf["hpa"] + canopy_load * 8000.0
    crown_ib = crown_ros_ftmin / 60.0 * crown_hpa * kw

    if o.engine_constants:       # D3, D4, D8: compare in Btu/ft/s and ft/min with the engine's constants
        i0_btu = I0 * KW_M_TO_BTU_FT_S_ENGINE
        tr = np.where(I0 > 1e-7, surf["ib"] / np.maximum(i0_btu, 1e-12), 0.0)
        cbd_eff = 16.0185 * np.asarray(cbd_kgm3, float) * kg2lb
        r0_ft = np.where(cbd_eff > 1e-7, FT_PER_M_ENGINE * 3.0 / np.maximum(cbd_eff, 1e-12), 0.0)
        ar = np.where(r0_ft > 1e-7, crown_ros_ftmin / np.maximum(r0_ft, 1e-12), 0.0)
    else:
        tr = np.where(I0 > 1e-7, I_s / I0, 0.0)
        ar = np.where(R0 > 1e-7, crown_ros / np.maximum(R0, 1e-7), 0.0)
    ftype = np.where(tr < 1.0, np.where(ar < 1.0, 0, 2), np.where(ar < 1.0, 1, 3))
    return dict(surface_ros=surf["ros"] / FT_PER_M, surface_I=I_s, I0=I0, R0=R0,
                crown_ros=crown_ros, crown_I=crown_ib, transition_ratio=tr,
                active_ratio=ar, fire_type=ftype, waf=waf)


def cruz2005_crown_ros(u10_kmh, cbd_kgm3, effm_pct):
    """
    Cruz, Alexander & Wakimoto (2005, Can. J. For. Res. 35:1626-1639) crown fire rate of spread.
      CROS_A = 11.02 * U10^0.90 * CBD^0.19 * exp(-0.17 * EFFM)        [m/min]   (eq. 3)
      CAC    = CROS_A / (3 / CBD)                                                 (eqs. 1-2)
      CROS_P = CROS_A * exp(-CAC)                                                 (eq. 4)
    U10: 10-m open wind (km/h); EFFM: estimated fine fuel moisture (%).
    Returns dict(active=CROS_A, passive=CROS_P, cac=CAC, ros=CROS_A if CAC>=1 else CROS_P).
    """
    u = np.maximum(np.asarray(u10_kmh, float), 1e-6)
    cbd = np.asarray(cbd_kgm3, float)
    a = 11.02 * u ** 0.90 * cbd ** 0.19 * np.exp(-0.17 * np.asarray(effm_pct, float))
    cac = a / (3.0 / cbd)
    p = a * np.exp(-cac)
    return dict(active=a, passive=p, cac=cac, ros=np.where(cac >= 1.0, a, p))
