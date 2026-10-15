// Behaviour of the engine at very steep slopes, wind = 0 (no wind term). With phi_s = 5.275 * beta^-0.3 * tan^2(slope),
// the engine's spread rate should satisfy ros = ros0 * (1 + phi_s), i.e. (ros/ros0 - 1) grows with tan^2(slope).
// SurfaceFire::calculateWindSpeedLimit() caps phi_s at 0.9 * reaction intensity (a wind speed in ft/min, numerically).
// Columns: slope %, tan^2, phi_s implied by the engine, phi_s predicted from the 100 % slope by tan^2 scaling, 0.9*IR.
#include <cstdio>
#include "behaveRun.h"
#include "fuelModels.h"
static double run(BehaveRun& b, int m, double slopePct, double* ir) {
  b.surface.updateSurfaceInputs(m, 6, 7, 8, 60, 90, FractionUnits::Percent, 0, SpeedUnits::FeetPerMinute,
    WindHeightInputMode::DirectMidflame, 0, WindAndSpreadOrientationMode::RelativeToUpslope, slopePct, SlopeUnits::Percent, 0,
    0, FractionUnits::Percent, 0, LengthUnits::Feet, 0, FractionUnits::Fraction);
  b.surface.doSurfaceRunInDirectionOfMaxSpread();
  if (ir) *ir = b.surface.getReactionIntensity(HeatSourceAndReactionIntensityUnits::BtusPerSquareFootPerMinute);
  return b.surface.getSpreadRate(SpeedUnits::FeetPerMinute);
}
int main() {
  FuelModels fm; SpeciesMasterTable sp; BehaveRun b(fm, sp);
  printf("fuel_model,slope_pct,tan2,phi_s_engine,phi_s_predicted_by_tan2,limit_0.9IR\n");
  for (int m : {1, 10, 165}) {
    double ir; double ros0 = run(b, m, 0, &ir); double phi100 = run(b, m, 100, nullptr) / ros0 - 1.0;
    for (double s : {100.0, 300.0, 1000.0, 3000.0, 10000.0, 30000.0}) {
      double t2 = (s / 100.0) * (s / 100.0); double phi = run(b, m, s, nullptr) / ros0 - 1.0;
      printf("%d,%g,%g,%.2f,%.2f,%.2f\n", m, s, t2, phi, phi100 * t2, 0.9 * ir);
    }
  }
}
