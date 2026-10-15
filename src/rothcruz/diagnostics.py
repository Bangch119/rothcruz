"""
Threshold-free diagnostics for comparing models with an observed spread rate over an interval.

required_crown_fraction  f = (R_obs - mean R_s) / (mean R_c - mean R_s): the fraction of the interval in active crown spread
                         needed to reproduce the observation; f > 1 means continuous crowning is insufficient.
wind_multiplier          the factor z on the wind for which a model's mean spread rate equals the observation,
                         e.g. the surface-only multiplier m_s.
"""
from __future__ import annotations

import math

import numpy as np


def required_crown_fraction(ros_obs, ros_surface, ros_crown):
    """Hourly (or per-record) surface and active crown spread rates are averaged over the interval first."""
    rs = float(np.mean(np.asarray(ros_surface, float)))
    rc = float(np.mean(np.asarray(ros_crown, float)))
    return (float(ros_obs) - rs) / (rc - rs) if rc != rs else math.nan


def wind_multiplier(ros_obs, mean_ros_at, lo=0.02, hi=1000.0, rtol=1e-10, max_iter=200):
    """Solve mean_ros_at(z) = ros_obs for z in [lo, hi] by bisection in log(z).

    mean_ros_at : function of the multiplier z returning the interval-mean spread rate; it must increase with z.
    Returns nan if there is no solution in the bracket (the model already exceeds the observation at `lo`, or stays
    below it at `hi`). The default bracket 0.02-1000 matches the evaluation in Bang and Bang."""
    g = lambda z: float(mean_ros_at(z)) - float(ros_obs)
    glo, ghi = g(lo), g(hi)
    if glo > 0 or ghi < 0:
        return math.nan
    a, b = math.log(lo), math.log(hi)
    for _ in range(max_iter):
        m = 0.5 * (a + b)
        if g(math.exp(m)) < 0: a = m
        else: b = m
        if b - a < rtol: break
    return math.exp(0.5 * (a + b))
