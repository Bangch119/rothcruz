#include <cstdio>
#include "behaveRun.h"
#include "fuelModels.h"
int main(){
  FuelModels fm; SpeciesMasterTable sp; BehaveRun b(fm, sp);
  int models[]={2,8,9,10,102,122,145,161,165,183,186,188,189};
  double m1s[]={3,6,10}; double w20s[]={0,5,10,20,30}; double cbhs[]={1,3,6}; double cbds[]={0.05,0.1,0.2}; double fmcs[]={80,100,120};
  double cc=50, ch=15, cr=0.5;
  printf("fm,m1,m10,m100,mlh,mlw,fmc,wind20_kmh,slope_pct,cbh_m,cbd,cc_pct,ch_m,cr,surf_ros_mmin,crown_ros_mmin,crown_ib_kwm,fire_type\n");
  for(int m: models) for(double m1: m1s) for(double w: w20s) for(double cbh: cbhs) for(double cbd: cbds) for(double fmc: fmcs){
    double m10=m1+1,m100=m1+2,mlh=60,mlw=90, slope=20;
    b.crown.updateCrownInputs(m,m1,m10,m100,mlh,mlw,fmc,FractionUnits::Percent,w,SpeedUnits::KilometersPerHour,
      WindHeightInputMode::TwentyFoot,0,WindAndSpreadOrientationMode::RelativeToUpslope,slope,SlopeUnits::Percent,0,
      cc,FractionUnits::Percent,ch,cbh,LengthUnits::Meters,cr,FractionUnits::Fraction,cbd,DensityUnits::KilogramsPerCubicMeter);
    b.crown.doCrownRunRothermel();
    printf("%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%.9g,%.9g,%.9g,%d\n",m,m1,m10,m100,mlh,mlw,fmc,w,slope,cbh,cbd,cc,ch,cr,
      b.crown.getSurfaceFireSpreadRate(SpeedUnits::MetersPerMinute), b.crown.getCrownFireSpreadRate(SpeedUnits::MetersPerMinute),
      b.crown.getCrownFirelineIntensity(FirelineIntensityUnits::KilowattsPerMeter), (int)b.crown.getFireType());
  }
}
