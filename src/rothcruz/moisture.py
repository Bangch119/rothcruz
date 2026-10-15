"""Dead fuel moisture from air temperature and relative humidity."""
from __future__ import annotations

import numpy as np


def simard_emc(rh_pct, temp_c):
    """Equilibrium moisture content (%) of Simard (1968), the regression on relative humidity (%) and air
    temperature (degrees F, converted here from degrees C) used in fire danger rating."""
    rh = np.asarray(rh_pct, float)
    tf = np.asarray(temp_c, float) * 9.0 / 5.0 + 32.0
    return np.where(rh < 10, 0.03229 + 0.281073 * rh - 0.000578 * rh * tf,
           np.where(rh < 50, 2.22749 + 0.160107 * rh - 0.01478 * tf,
                    21.0606 + 0.005565 * rh ** 2 - 0.00035 * rh * tf - 0.483199 * rh))


def dead_fuel_moisture(rh_pct, temp_c, factor=1.03, step=(1.0, 2.0)):
    """1-, 10- and 100-h dead fuel moisture (%): factor * EMC for 1-h fuels, plus `step` percentage points for the
    10-h and 100-h classes (the convention of the evaluation in Bang and Bang, after Rothermel 1991 and Cruz et al. 2005)."""
    m1 = factor * simard_emc(rh_pct, temp_c)
    return m1, m1 + step[0], m1 + step[1]
