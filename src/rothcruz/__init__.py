"""rothcruz: Rothermel surface fire and crown fire spread models verified against the BehavePlus engine."""
from .options import BEHAVE_COMPAT, PUBLISHED, ModelOptions
from .surface import FUEL_MODELS, fuel_model, rothermel, surface_fire
from .crown import (crown_fire, cruz2005_crown_ros, critical_intensity, critical_active_ros, wind_adjustment_factor)
from .scott_reinhardt import scott_reinhardt
from .moisture import dead_fuel_moisture, simard_emc
from .diagnostics import required_crown_fraction, wind_multiplier

__version__ = "0.1.0"
__all__ = ["ModelOptions", "PUBLISHED", "BEHAVE_COMPAT", "FUEL_MODELS", "fuel_model", "rothermel", "surface_fire", "crown_fire", "cruz2005_crown_ros", "critical_intensity",
           "critical_active_ros", "wind_adjustment_factor", "scott_reinhardt", "simard_emc", "dead_fuel_moisture",
           "required_crown_fraction", "wind_multiplier", "__version__"]
