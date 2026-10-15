# rothcruz

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23387346.svg)](https://doi.org/10.5281/zenodo.23387346)

Rothermel surface fire spread model and crown fire models in Python, verified case by case against the
official BehavePlus calculation engine ([firelab/behave](https://github.com/firelab/behave)).

**Version 0.2.0.** The verification suite passes and the `published` / `behave_compat` modes work. Documentation is
this README and `docs/divergences.md`.

## What it computes

* Surface fire: Rothermel (1972) with the revisions of Albini (1976), the 53 standard fuel models, dynamic curing
  (Scott and Burgan 2005), wind adjustment factor.
* Crown fire: Van Wagner (1977) initiation and active-crown criteria; crown fire spread rate of Rothermel (1991) and of
  Cruz et al. (2005) (active and passive); crown fraction burned and final spread rate of Scott and Reinhardt (2001).
* Dead fuel moisture from temperature and relative humidity (equilibrium moisture content of Simard 1968).
* Diagnostics for comparing models with an observed spread rate: the required fraction of time in active crown spread
  `required_crown_fraction` and the wind multiplier `wind_multiplier` (e.g. the surface-only multiplier m_s).
* All functions take scalars or NumPy arrays.

`python examples/uiseong_2025_interval.py` reproduces the worked example of the software paper (one interval of the
2025 Uiseong fire, Korea).

## Install and use

```
pip install -e .            # from a clone of this repository
```

```python
import numpy as np
from rothcruz import BEHAVE_COMPAT, crown_fire, cruz2005_crown_ros

wind10 = np.array([20.0, 30.0, 44.0])            # 10-m open wind (km/h)
site = dict(fuel=165,                            # fuel model TU5
            m1=4, m10=5, m100=6, mlh=30, mlw=90, # fuel moisture (%)
            fmc=100, slope_pct=20.0,             # foliar moisture, slope (%)
            cbh_m=2.0, cbd_kgm3=0.15,            # canopy base height, CBD
            canopy_cover_pct=50, canopy_height_m=15.0)
u20 = wind10 / 1.15                              # 20-ft wind (km/h)
out = crown_fire(wind20_kmh=u20, **site)         # published mode (default)
r91 = out["crown_ros"]                           # Rothermel (1991), m/min
c05 = cruz2005_crown_ros(wind10, 0.15, 4.0)["active"]  # Cruz et al. (2005)
eng = crown_fire(wind20_kmh=u20, **site, options=BEHAVE_COMPAT)  # as the engine
```


## Two modes

Every function takes `options=`. `PUBLISHED` (the default) follows the original equations and tables; `BEHAVE_COMPAT`
reproduces the BehavePlus engine. Single differences can be switched with `ModelOptions`; each field is one entry of
[docs/divergences.md](docs/divergences.md).

```python
from rothcruz import BEHAVE_COMPAT, PUBLISHED, crown_fire

out_engine = crown_fire(165, 4, 5, 6, 30, 90, 100, 25.0, 20.0, 2.0, 0.15, 50, 15.0, options=BEHAVE_COMPAT)
out_limit  = crown_fire(165, 4, 5, 6, 30, 90, 100, 25.0, 20.0, 2.0, 0.15, 50, 15.0,
                        options=PUBLISHED.replace(wind_limit="effective"))   # Rothermel (1972) wind limit
```

| Field | Engine behaviour reproduced when set | ID |
|---|---|---|
| `curing="behaveplus"` | rounded curing coefficients | D1 |
| `override_crown_ratio=True` | a user crown ratio is replaced by (CH - CBH)/CH | D2 |
| `engine_constants=True` | the engine's unit constants and classification order | D3, D4, D8 |
| `van_wagner_limits=True` | FMC at least 30 %, CBH at least 0.1 m | D5 |
| `wind_limit` | `"none"` (both presets), `"midflame"`, `"effective"` | D6 |
| `cap_slope_factor=True` | slope factor capped at 0.9 I_R (very steep slopes) | D7 |
| `fuel_table="engine"` | live SAV of fuel model 7 is 1500, not Albini's 1550 | D9 |
| `crown_slope=False` | crown fire spread rate without slope (also Rothermel's rule for runs across uneven terrain) | D10 |

## Verification against the engine

`pip install -e ".[test]" && pytest` compares the code with outputs of the engine at a pinned commit
(`tests/reference/README.md` explains how they were generated). The reference data, drivers and a script that regenerates
them byte for byte are archived as a dataset: https://doi.org/10.5281/zenodo.23386397

| Quantity | Cases | Largest relative difference |
|---|---|---|
| Surface fire spread rate | 7,632 | 2.6e-6 |
| Crown fire spread rate (Rothermel 1991) | 5,265 | 6.1e-9 |
| Scott and Reinhardt final spread rate | 33,696 | 4.5e-6 |
| Fire type (surface, torching, conditional, crowning) | 38,961 | identical in all cases |
| Slope factor cap at very steep slopes | 18 | within 1e-3 (the reference prints two decimals) |

All of these use `BEHAVE_COMPAT`. The wind-limit options cannot be compared with the engine (it cannot enable the limit);
they are tested analytically (`tests/test_modes.py`).

## AI assistance

Parts of the code, tests and documentation were developed with the assistance of an AI model.
The authors reviewed the code, verified it against the engine outputs above, and are responsible for it.

## Citation and license

Cite the archived release: https://doi.org/10.5281/zenodo.23387346 (metadata in `CITATION.cff`). BSD 3-Clause, see `LICENSE.txt`.
