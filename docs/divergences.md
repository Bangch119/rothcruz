# Where the BehavePlus engine differs from the original publications

Engine: firelab/behave commit `d4acc8a` (2026-09-14). Status = the field of `ModelOptions` that controls it. `PUBLISHED` (default) follows the publications; `BEHAVE_COMPAT` reproduces the engine.

| ID | Item | Published | Engine | Effect | Status here |
|---|---|---|---|---|---|
| D1 | Dynamic curing fraction | (1.20 − M_lh)/0.90 (Scott & Burgan 2005) | 1.333 − 1.11 M_lh (rounded; exact form kept as a comment, `surfaceFuelbedIntermediates.cpp` 646–647) | up to 6.5 % in fireline intensity on the verification grid, about 12 % between grid points | `ModelOptions.curing`: `"exact"` (published) or `"behaveplus"` |
| D2 | Crown ratio in crown fire runs | user input | replaced by (stand height − CBH)/stand height when the wind adjustment factor uses the canopy (`crown.cpp` 112–122, 172–182) | changes midflame wind and surface spread rate | `override_crown_ratio` (True in `BEHAVE_COMPAT`) |
| D3 | Btu to kJ | 1.055056 kJ (3.46147 kW/m per Btu/ft/s) | 1.05587 kJ (3.464140419) in most conversions; reaction intensity uses the IT value (`behaveUnits.cpp`) | kW/m outputs 0.08 % high | `engine_constants` (True in `BEHAVE_COMPAT`) |
| D4 | Crown transition ratio | I_s compared with I0 in kW/m | I0 converted to Btu/ft/s with 0.2886719, then compared | crown-initiation boundary moves by 0.077 % | `engine_constants` (True in `BEHAVE_COMPAT`) |
| D5 | Van Wagner limits | none | FMC ≥ 30 %, CBH ≥ 0.1 m (`crown.cpp` 493–509) | changes I0 when below the limits | `van_wagner_limits` (True in `BEHAVE_COMPAT`) |
| D6 | Wind limit | effective wind ≤ 0.9 I_R | off by default since e6d2b8a (2025-07-21) and cannot be enabled: flag reset at every calculation, setter protected and unused (firelab/behave issue 77) | spread rate unlimited at high wind | `wind_limit`: `"none"` (default in both presets), `"midflame"` or `"effective"` (Rothermel 1972); tested analytically, the engine cannot enable it |
| D7 | Slope factor cap | none | φ_s capped at 0.9 I_R (a wind speed) for very steep slopes, about 77° and above (`surfaceFire.cpp` 465–475) | none on normal slopes | `cap_slope_factor` (True in `BEHAVE_COMPAT`), tested against the engine up to 30,000 % slope |
| D8 | Density conversion | 0.0624279606 lb/ft³ per kg/m³ | 0.06242781786 (1/16.0185) | 2.3e-6 relative; moves O'active by 1.7e-6 and the crown fraction burned by up to 2e-5 | `engine_constants` (True in `BEHAVE_COMPAT`) |
| D9 | Fuel model 7, live SAV | 1550 ft⁻¹ (Albini 1976, Table 7, p. 92) | 1500 ft⁻¹ (firelab/behave issue 76) | spread rate of FM7 up to 1.5 % (median 0.7 %) | `fuel_table`: `"published"` (1550, default) or `"engine"` (1500) |
| D10 | Slope in the crown fire spread rate | slope factor of fuel model 10 included (Rothermel 1991, Table 1a and nomograms; Scott & Reinhardt 2001, Eq. 8, 19, 20). Rothermel (1991, p. 25) recommends zero slope for the average rate of a run across uneven terrain | zero slope in the fuel model 10 run and in O'active | crown fire rate higher on slopes: +3.5 to +64 % (median +17 %) over the 19 intervals of Bang and Bang; none on flat ground | `crown_slope`: True in PUBLISHED, False in BEHAVE_COMPAT; O'active tested for consistency with R'active |

Not a divergence: the engine hard-codes the propagating flux ratio and the wind factors B and K of fuel model 10 used for
O'active; they equal the computed values to 1e-16.
