"""
Model options: one object that says, for every documented divergence between the original publications and the
BehavePlus engine, which behaviour is used. Two presets cover the common cases:

    PUBLISHED      the original equations and tables (default)
    BEHAVE_COMPAT  reproduces the engine (firelab/behave commit d4acc8a)

Each field corresponds to an entry of docs/divergences.md (D1-D10).
"""
from __future__ import annotations

from dataclasses import dataclass, replace

CURING = ("exact", "behaveplus")
WIND_LIMITS = ("none", "midflame", "effective")
FUEL_TABLES = ("published", "engine")


@dataclass(frozen=True)
class ModelOptions:
    curing: str = "exact"                  # D1  "exact": (1.20 - M)/0.90   "behaveplus": 1.333 - 1.11 M
    wind_limit: str = "none"               # D6  "none" | "midflame" (cap midflame wind) | "effective" (cap effective wind, Rothermel 1972)
    fuel_table: str = "published"          # D9  "published": Albini (1976) live SAV of FM7 = 1550   "engine": 1500
    override_crown_ratio: bool = False     # D2  True: ignore a user crown ratio and use (CH - CBH)/CH, as the engine does
    van_wagner_limits: bool = False        # D5  True: FMC >= 30 %, CBH >= 0.1 m in the critical intensity, as the engine does
    crown_slope: bool = True               # D10 True: slope factor of fuel model 10 in the crown fire spread rate and O'active
                                           #     (Rothermel 1991, Table 1a; Scott & Reinhardt 2001, Eq. 8, 19); the engine uses zero slope.
                                           #     Rothermel (1991, p. 25) recommends zero slope for the average rate of runs across uneven terrain.
    cap_slope_factor: bool = False         # D7  True: cap phi_s at 0.9 * reaction intensity (engine), matters above about 77 degrees
    engine_constants: bool = False         # D3, D4, D8  True: the engine's rounded unit constants and its classification order

    def __post_init__(self):
        if self.curing not in CURING:
            raise ValueError(f"curing must be one of {CURING}, got {self.curing!r}")
        if self.wind_limit not in WIND_LIMITS:
            raise ValueError(f"wind_limit must be one of {WIND_LIMITS}, got {self.wind_limit!r}")
        if self.fuel_table not in FUEL_TABLES:
            raise ValueError(f"fuel_table must be one of {FUEL_TABLES}, got {self.fuel_table!r}")

    def replace(self, **changes) -> "ModelOptions":
        """Copy with some fields changed, e.g. PUBLISHED.replace(wind_limit="effective")."""
        return replace(self, **changes)


PUBLISHED = ModelOptions()
BEHAVE_COMPAT = ModelOptions(curing="behaveplus", wind_limit="none", fuel_table="engine", override_crown_ratio=True,
                             van_wagner_limits=True, cap_slope_factor=True, engine_constants=True, crown_slope=False)


def resolve(options) -> ModelOptions:
    """None -> PUBLISHED; "published" / "behave_compat" -> the preset; a ModelOptions is returned unchanged."""
    if options is None:
        return PUBLISHED
    if isinstance(options, ModelOptions):
        return options
    if isinstance(options, str):
        presets = {"published": PUBLISHED, "behave_compat": BEHAVE_COMPAT}
        if options in presets:
            return presets[options]
        raise ValueError(f"unknown mode {options!r}; use 'published', 'behave_compat' or a ModelOptions")
    raise TypeError("options must be None, a mode name or a ModelOptions")
