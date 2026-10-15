"""
Worked example (Section 3.1 of the software paper): the run of the 2025 Uiseong fire, Korea, from 14:29 on 25 March to
02:49 on 26 March (KST), mean head-fire spread rate 66.06 m/min from VIIRS active-fire detections (Bang and Bang).

Weather: hourly records of the Korea Meteorological Administration (ASOS, Open MET Data Portal) at Yeongdeok and
Cheongsong. Cheongsong has no record from 21:00 to 02:00; as in the evaluation, the first three hours of that gap are
filled by linear interpolation (values with more decimals below) and the last three are left out, giving 21 station-hours.
Run:  python examples/uiseong_2025_interval.py
"""
import numpy as np
from rothcruz import (PUBLISHED, crown_fire, cruz2005_crown_ros, dead_fuel_moisture, required_crown_fraction,
                      scott_reinhardt, wind_multiplier)

ROS_OBS = 66.06465800448147        # m/min, observed (VIIRS)
SLOPE = 18.763040091685248         # %, slope along the run direction (root mean square of upslope components)
HOURS = [  # station, KST, 10-m wind (m/s), air temperature (C), relative humidity (%)
    ("Yeongdeok", "2025-03-25 15:00", 4.7, 25.8, 17),
    ("Yeongdeok", "2025-03-25 16:00", 4, 25.7, 18),
    ("Yeongdeok", "2025-03-25 17:00", 2.6, 24.9, 19),
    ("Yeongdeok", "2025-03-25 18:00", 3.4, 25.7, 19),
    ("Yeongdeok", "2025-03-25 19:00", 4.8, 25, 21),
    ("Yeongdeok", "2025-03-25 20:00", 8.8, 21, 32),
    ("Yeongdeok", "2025-03-25 21:00", 10.5, 18.7, 29),
    ("Yeongdeok", "2025-03-25 22:00", 12.2, 17.1, 26),
    ("Yeongdeok", "2025-03-25 23:00", 10.5, 16.3, 25),
    ("Yeongdeok", "2025-03-26 00:00", 4.1, 15.1, 28),
    ("Yeongdeok", "2025-03-26 01:00", 2.2, 13.6, 31),
    ("Yeongdeok", "2025-03-26 02:00", 2.6, 13.3, 33),
    ("Cheongsong", "2025-03-25 15:00", 3.1, 25.1, 16),
    ("Cheongsong", "2025-03-25 16:00", 4.6, 26, 12),
    ("Cheongsong", "2025-03-25 17:00", 5.4, 25.6, 13),
    ("Cheongsong", "2025-03-25 18:00", 8.7, 28, 15),
    ("Cheongsong", "2025-03-25 19:00", 8.5, 22.9, 23),
    ("Cheongsong", "2025-03-25 20:00", 7, 19.4, 29),
    ("Cheongsong", "2025-03-25 21:00", 6.28571, 18.1143, 30.1429),
    ("Cheongsong", "2025-03-25 22:00", 5.57143, 16.8286, 31.2857),
    ("Cheongsong", "2025-03-25 23:00", 4.85714, 15.5429, 32.4286),
]
SITE = dict(fuel=165, mlh=30, mlw=90, fmc=100, slope_pct=SLOPE, cbh_m=2.0, cbd_kgm3=0.15,
            canopy_cover_pct=50, canopy_height_m=15.0)   # TU5 and the base canopy of the evaluation


# The run crossed ridges and valleys, so the crown fire spread rate uses zero slope, following Rothermel (1991, p. 25):
# "To estimate the average rate of spread wherein the fire may travel across uneven terrain, use zero slope."
# The surface fire still uses the slope of the run. This is how the evaluation of Bang and Bang was run.
OPTIONS = PUBLISHED.replace(crown_slope=False)


def main(options=OPTIONS, verbose=True):
    wind_ms = np.array([h[2] for h in HOURS]); temp = np.array([h[3] for h in HOURS]); rh = np.array([h[4] for h in HOURS])
    m1, m10, m100 = dead_fuel_moisture(rh, temp)
    u20 = wind_ms * 3.6 / 1.15                      # 10-m wind (m/s) -> 20-ft wind (km/h)
    cf = crown_fire(m1=m1, m10=m10, m100=m100, wind20_kmh=u20, options=options, **SITE)
    c05 = cruz2005_crown_ros(wind_ms * 3.6, 0.15, m1)["active"]
    sr = scott_reinhardt(m1=m1, m10=m10, m100=m100, wind20_kmh=u20, options=options, **SITE)
    surface_at = lambda z: np.mean(crown_fire(m1=m1, m10=m10, m100=m100, wind20_kmh=z * u20, options=options, **SITE)["surface_ros"])
    out = dict(n=len(HOURS), surface=float(np.mean(cf["surface_ros"])), r91=float(np.mean(cf["crown_ros"])),
               c05=float(np.mean(c05)), scott_reinhardt=float(np.mean(sr["final_ros"])), cfb=float(np.mean(sr["cfb"])),
               f_r91=required_crown_fraction(ROS_OBS, cf["surface_ros"], cf["crown_ros"]),
               f_c05=required_crown_fraction(ROS_OBS, cf["surface_ros"], c05),
               m_s=wind_multiplier(ROS_OBS, surface_at),
               torching=int(np.sum(sr["fire_type"] == 1)), crowning=int(np.sum(sr["fire_type"] == 3)))
    if verbose:
        print(f"{out['n']} station-hours; observed {ROS_OBS:.2f} m/min")
        print(f"mean spread rate (m/min): surface {out['surface']:.2f}, Rothermel (1991) crown {out['r91']:.2f}, "
              f"Cruz et al. (2005) crown {out['c05']:.2f}, Scott & Reinhardt final {out['scott_reinhardt']:.2f}")
        print(f"required crown fraction f: {out['f_r91']:.2f} (Rothermel 1991), {out['f_c05']:.2f} (Cruz et al. 2005)")
        print(f"surface-only wind multiplier m_s = {out['m_s']:.1f}; mean crown fraction burned {out['cfb']:.2f}; "
              f"torching {out['torching']}, crowning {out['crowning']} of {out['n']} station-hours")
    return out


if __name__ == "__main__":
    main()
