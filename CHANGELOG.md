# Changelog

## 0.2.0
- New option `crown_slope` (D10). Rothermel (1991, Table 1a) and Scott & Reinhardt (2001, Eq. 8, 19, 20) include the slope
  factor of fuel model 10 in the crown fire spread rate and in O'active; the BehavePlus engine uses zero slope.
  `PUBLISHED` now includes the slope (`crown_slope=True`); `BEHAVE_COMPAT` keeps zero slope. Results on sloped ground change
  in published mode. Rothermel (1991, p. 25) recommends zero slope for the average rate of a run across uneven terrain;
  use `PUBLISHED.replace(crown_slope=False)` for that, as the worked example does.
- O'active is 0 (instead of undefined) when the slope alone makes the crown fire rate reach R'active.

## 0.1.0
- First release: surface fire, crown fire (Van Wagner, Rothermel 1991, Cruz et al. 2005, Scott & Reinhardt 2001),
  PUBLISHED and BEHAVE_COMPAT modes (D1-D9), dead fuel moisture, diagnostics, verification against the engine.
