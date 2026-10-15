"""
Rothermel (1972) surface fire spread model with Albini (1976) revisions and
Scott & Burgan (2005) dynamic fuel-model curing.

Implemented independently from the published equations, vectorised with NumPy
so that every input may be a scalar or an array (raster maps).

References
----------
Rothermel, R.C. (1972) USDA FS Res. Pap. INT-115.
Albini, F.A. (1976) USDA FS Gen. Tech. Rep. INT-30.
Scott, J.H. & Burgan, R.E. (2005) USDA FS Gen. Tech. Rep. RMRS-GTR-153.
Byram, G.M. (1959) Combustion of forest fuels. In: Forest Fire: Control and Use.

Internal units follow the original papers (English units: ft, lb, Btu, min).
The public function `surface_fire` accepts and returns SI units.
"""
from __future__ import annotations

import csv
import os
import dataclasses
from importlib.resources import files
from dataclasses import dataclass

import numpy as np

from .constants import BTU_FT_S_TO_KW_M_ENGINE
from .options import resolve

# ----------------------------------------------------------------- constants
RHO_P = 32.0        # oven-dry particle density, lb/ft3
S_T = 0.0555        # total mineral content
S_E = 0.010         # effective (silica-free) mineral content
SAV_10H = 109.0     # 1/ft
SAV_100H = 30.0     # 1/ft

FT_PER_M = 3.280839895
MPH_TO_FTMIN = 88.0
BTU_FT_S_TO_KW_M = 1.05505585262 / 0.3048  # = 3.46147 kW/m (IT Btu). NB: BehavePlus uses 3.464140
BTU_FT2_TO_KJ_M2 = 11.3565   # 1 Btu/ft2 = 11.3565 kJ/m2
BTU_FT2MIN_TO_KW_M2 = 0.189275

_DATA = str(files("rothcruz") / "data" / "fuel_models.csv")


@dataclass(frozen=True)
class FuelModel:
    number: int
    code: str
    depth: float          # ft
    mx_dead: float        # fraction
    h_dead: float         # Btu/lb
    h_live: float         # Btu/lb
    w: tuple              # (1h, 10h, 100h, live herb, live woody) lb/ft2
    sav1: float
    savherb: float
    savwoody: float
    dynamic: bool


def load_fuel_models(path: str = _DATA) -> dict:
    """Return {number: FuelModel, code: FuelModel} for the 53 standard models."""
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            fm = FuelModel(
                number=int(r["number"]), code=r["code"],
                depth=float(r["depth_ft"]), mx_dead=float(r["mx_dead"]),
                h_dead=float(r["h_dead"]), h_live=float(r["h_live"]),
                w=tuple(float(r[k]) for k in ("w1", "w10", "w100", "wherb", "wwoody")),
                sav1=float(r["sav1"]), savherb=float(r["savherb"]),
                savwoody=float(r["savwoody"]), dynamic=bool(int(r["dynamic"])),
            )
            out[fm.number] = fm
            out[fm.code] = fm
    return out


FUEL_MODELS = load_fuel_models()      # tables as in the BehavePlus engine (extracted from its source)

# Entries where the published tables differ from the engine's (D9)
_PUBLISHED_OVERRIDES = {7: dict(savwoody=1550.0)}      # Albini (1976), Table 7, p. 92: live SAV of model 7 is 1550 ft-1


def fuel_model(key, options=None) -> "FuelModel":
    """Fuel model by number or code, with the table that `options.fuel_table` selects ("published" or "engine")."""
    o = resolve(options)
    if isinstance(key, FuelModel):                      # a model supplied by the user is used as it is
        return key
    fm = FUEL_MODELS[key]
    if o.fuel_table == "published" and fm.number in _PUBLISHED_OVERRIDES:
        return dataclasses.replace(fm, **_PUBLISHED_OVERRIDES[fm.number])
    return fm


def _size_class(sav):
    """Albini (1976) FIREMOD size-class bins used for the g_ij weights."""
    return np.select([sav >= 1200, sav >= 192, sav >= 96, sav >= 48, sav >= 16],
                     [1, 2, 3, 4, 5], 6)


def _safe_div(a, b):
    a, b = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float))
    out = np.zeros_like(a)
    m = b > 0
    out[m] = a[m] / b[m]
    return out


def rothermel(fm: FuelModel, m1, m10, m100, mlh, mlw,
              wind_midflame_ftmin, slope_tan, *, options=None):
    """
    Core model in original units. Moistures are fractions (0.06 = 6 %).
    Wind is midflame speed in ft/min, assumed aligned with the upslope
    direction (head-fire, direction of maximum spread).

    options : ModelOptions, a mode name ("published", "behave_compat") or None (= published).
              Used here: curing ("exact" (1.20 - Mlh)/0.90 or "behaveplus" 1.333 - 1.11*Mlh),
              wind_limit ("none", "midflame": cap the midflame wind at 0.9*I_R, "effective": cap the effective wind,
              as in Rothermel 1972) and cap_slope_factor.  The fuel model itself is chosen by the caller
              (see fuel_model()).

    Returns dict with ros (ft/min), ir (Btu/ft2/min), hpa (Btu/ft2),
    ib (Btu/ft/s), flame (ft), plus intermediates.
    """
    o = resolve(options)
    m1, m10, m100, mlh, mlw, U, tanS = np.broadcast_arrays(
        *[np.asarray(x, dtype=float) for x in
          (m1, m10, m100, mlh, mlw, wind_midflame_ftmin, slope_tan)])
    shape = m1.shape

    # --- particle classes: dead [1h,10h,100h,dead-herb], live [herb, woody]
    w1, w10, w100, wh, ww = fm.w
    if fm.dynamic and wh > 0:
        # Scott & Burgan (2005): fully cured at Mlh <= 30 %, uncured at >= 120 %
        if o.curing == "behaveplus":
            cured = np.where(mlh < 0.30, 1.0,
                             np.where(mlh <= 1.20, 1.333 - 1.11 * mlh, 0.0))
        else:
            cured = np.clip((1.20 - mlh) / 0.90, 0.0, 1.0)
    else:
        cured = np.zeros(shape)
    w_d = [np.full(shape, w1), np.full(shape, w10), np.full(shape, w100), wh * cured]
    w_l = [wh * (1.0 - cured), np.full(shape, ww)]
    sav_d = [fm.sav1, SAV_10H, SAV_100H, fm.savherb]
    sav_l = [fm.savherb, fm.savwoody]
    mf_d = [m1, m10, m100, m1]          # dead herb takes 1-h moisture
    mf_l = [mlh, mlw]

    # --- surface-area weights (Rothermel 1972 eq. 53-57)
    a_d = [s * w / RHO_P for s, w in zip(sav_d, w_d)]
    a_l = [s * w / RHO_P for s, w in zip(sav_l, w_l)]
    A_d, A_l = sum(a_d), sum(a_l)
    A_T = A_d + A_l
    f_d = [_safe_div(a, A_d) for a in a_d]
    f_l = [_safe_div(a, A_l) for a in a_l]
    F_d, F_l = _safe_div(A_d, A_T), _safe_div(A_l, A_T)

    # --- g_ij: weights summed over Albini size classes (for net loading)
    def g_weights(f, sav):
        cls = [_size_class(np.float64(s)) if s > 0 else 0 for s in sav]
        return [sum(f[k] for k in range(len(f)) if cls[k] == cls[j] and sav[k] > 0)
                if sav[j] > 0 else np.zeros(shape) for j in range(len(f))]
    g_d, g_l = g_weights(f_d, sav_d), g_weights(f_l, sav_l)

    # --- characteristic SAV, packing ratio
    sig_d = sum(f * s for f, s in zip(f_d, sav_d))
    sig_l = sum(f * s for f, s in zip(f_l, sav_l))
    sigma = F_d * sig_d + F_l * sig_l
    w_tot = sum(w_d) + sum(w_l)
    beta = w_tot / RHO_P / fm.depth
    beta_op = np.where(sigma > 0, 3.348 * np.power(np.maximum(sigma, 1e-9), -0.8189), 1.0)
    rel = _safe_div(beta, beta_op)

    # --- net loading, heat, moisture per life category
    wn_d = sum(g * w * (1 - S_T) for g, w in zip(g_d, w_d))
    wn_l = sum(g * w * (1 - S_T) for g, w in zip(g_l, w_l))
    Mf_d = sum(f * m for f, m in zip(f_d, mf_d))
    Mf_l = sum(f * m for f, m in zip(f_l, mf_l))

    # live moisture of extinction (Albini 1976 revision of eq. 88)
    eps_d = [w * np.exp(-138.0 / s) if s > 0 else 0 * w for w, s in zip(w_d, sav_d)]
    eps_l = [w * np.exp(-500.0 / s) if s > 0 else 0 * w for w, s in zip(w_l, sav_l)]
    Wd, Wl = sum(eps_d), sum(eps_l)
    mfdead_fine = _safe_div(sum(e * m for e, m in zip(eps_d, mf_d)), Wd)
    mx_live = np.where((Wd > 0) & (Wl > 0),
                       2.9 * _safe_div(Wd, Wl) * (1 - mfdead_fine / fm.mx_dead) - 0.226,
                       fm.mx_dead)
    mx_live = np.maximum(mx_live, fm.mx_dead)

    def eta_m(mf, mx):
        r = np.clip(_safe_div(mf, mx), 0, 1)
        r = np.where(mx > 0, r, 1.0)
        return 1 - 2.59 * r + 5.11 * r ** 2 - 3.52 * r ** 3
    eta_s = min(0.174 * S_E ** -0.19, 1.0)

    # --- reaction velocity (Albini A-exponent)
    s15 = np.power(sigma, 1.5)
    gmax = s15 / (495.0 + 0.0594 * s15)
    A = np.where(sigma > 0, 133.0 * np.power(np.maximum(sigma, 1e-9), -0.7913), 0.0)
    gamma = gmax * np.power(rel, A) * np.exp(A * (1.0 - rel))

    ir = gamma * (wn_d * fm.h_dead * eta_m(Mf_d, fm.mx_dead) * eta_s * (A_d > 0)
                  + wn_l * fm.h_live * eta_m(Mf_l, mx_live) * eta_s * (A_l > 0))

    xi = np.exp((0.792 + 0.681 * np.sqrt(sigma)) * (beta + 0.1)) / (192.0 + 0.2595 * sigma)

    # --- wind & slope factors
    C = 7.47 * np.exp(-0.133 * np.power(sigma, 0.55))
    B = 0.02526 * np.power(sigma, 0.54)
    E = 0.715 * np.exp(-3.59e-4 * sigma)
    u_lim = 0.9 * ir                                   # Rothermel (1972) wind speed limit (ft/min)
    if o.wind_limit == "midflame":
        U = np.minimum(U, u_lim)
    rel_term = np.power(np.maximum(rel, 1e-12), -E)
    phi_w = C * np.power(U, B) * rel_term
    phi_s = 5.275 * np.power(np.maximum(beta, 1e-12), -0.3) * tanS ** 2
    if o.cap_slope_factor:                              # D7: the engine compares phi_s with the wind limit
        phi_s = np.where(phi_s > 0, np.minimum(phi_s, u_lim), phi_s)
    phi_e = phi_w + phi_s
    wind_limited = np.zeros(shape, bool)
    if o.wind_limit == "effective":                     # limit the effective wind speed: phi_E <= C * u_lim^B * rel^-E
        cap = C * np.power(np.maximum(u_lim, 0.0), B) * rel_term
        wind_limited = phi_e > cap
        phi_e = np.minimum(phi_e, cap)

    # --- heat sink
    rho_b = w_tot / fm.depth
    def qsum(f, sav, mf):
        return sum(fi * (np.exp(-138.0 / s) if s > 0 else 0.0) * (250 + 1116 * m) * (m > 0)
                   for fi, s, m in zip(f, sav, mf))
    sink = rho_b * (F_d * qsum(f_d, sav_d, mf_d) + F_l * qsum(f_l, sav_l, mf_l))

    ros = _safe_div(ir * xi * (1 + phi_e), sink)
    t_res = _safe_div(384.0, sigma)                # min
    hpa = ir * t_res                               # Btu/ft2
    ib = hpa * ros / 60.0                          # Btu/ft/s
    flame = 0.45 * np.power(ib, 0.46)              # ft (Byram 1959)
    return dict(ros=ros, ir=ir, hpa=hpa, ib=ib, flame=flame,
                sigma=sigma, beta=beta, beta_op=beta_op, xi=xi,
                phi_w=phi_w, phi_s=phi_s, phi_e=phi_e, wind_limited=wind_limited, mx_live=mx_live,
                sink=sink, wind_C=C, wind_B=B, wind_E=E, rel_packing=rel)


def surface_fire(fuel, m1, m10, m100, mlh, mlw, wind_midflame_ms, slope_deg=0.0, *, options=None):
    """
    SI wrapper. Moisture in %, midflame wind in m/s, slope in degrees.
    Returns ros (m/min), I_B (kW/m), flame length (m), HPA (kJ/m2), IR (kW/m2).
    """
    o = resolve(options)
    fm = fuel_model(fuel, o)
    kw = BTU_FT_S_TO_KW_M_ENGINE if o.engine_constants else BTU_FT_S_TO_KW_M
    pct = lambda x: np.asarray(x, float) / 100.0
    U = np.asarray(wind_midflame_ms, float) * FT_PER_M * 60.0
    tanS = np.tan(np.radians(np.asarray(slope_deg, float)))
    r = rothermel(fm, pct(m1), pct(m10), pct(m100), pct(mlh), pct(mlw), U, tanS, options=o)
    return dict(ros=r["ros"] / FT_PER_M, I_B=r["ib"] * kw,
                flame=r["flame"] / FT_PER_M, hpa=r["hpa"] * BTU_FT2_TO_KJ_M2,
                ir=r["ir"] * BTU_FT2MIN_TO_KW_M2, _raw=r)
