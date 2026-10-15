"""Unit constants of the BehavePlus engine that differ slightly from the exact values (see docs/divergences.md)."""
KW_M_TO_BTU_FT_S_ENGINE = 0.2886719      # I0 (kW/m) -> Btu/ft/s before the transition test (D4)
BTU_FT_S_TO_KW_M_ENGINE = 3.464140419    # 1.05587 kJ per Btu in most conversions; the exact (IT) value is 3.46147 (D3)
KGM3_TO_LBFT3_ENGINE = 0.06242781786     # 1 / 16.0185; the exact value is 0.0624279606 (D8)
FT_PER_M_ENGINE = 3.28084
