#include <iostream>
#include <cstdio>
#include "behaveRun.h"
#include "fuelModels.h"
int main(){
  FuelModels fm; SpeciesMasterTable sp; BehaveRun b(fm, sp);
  int models[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,
    101,102,103,104,105,106,107,108,109,121,122,123,124,
    141,142,143,144,145,146,147,148,149,161,162,163,164,165,
    181,182,183,184,185,186,187,188,189,201,202,203,204};
  double m1s[]={3,6,10}; double herbs[]={30,60,90,120}; double winds[]={0,2,5,10}; double slopes[]={0,30,60};
  printf("fm,m1,m10,m100,mlh,mlw,wind_mph,slope_pct,ros_ftmin,ib_btufts,fl_ft,ir_btuft2min\n");
  for(int m: models){ if(!fm.isFuelModelDefined(m)) {fprintf(stderr,"undef %d\n",m); continue;}
   for(double m1: m1s) for(double mh: herbs) for(double w: winds) for(double s: slopes){
    double m10=m1+1, m100=m1+2, mlw=mh+30;
    b.surface.updateSurfaceInputs(m,m1,m10,m100,mh,mlw,FractionUnits::Percent,w,SpeedUnits::MilesPerHour,
      WindHeightInputMode::DirectMidflame,0,WindAndSpreadOrientationMode::RelativeToUpslope,s,SlopeUnits::Percent,0,
      0,FractionUnits::Percent,0,LengthUnits::Feet,0,FractionUnits::Fraction);
    b.surface.doSurfaceRunInDirectionOfMaxSpread();
    printf("%d,%g,%g,%g,%g,%g,%g,%g,%.9g,%.9g,%.9g,%.9g\n",m,m1,m10,m100,mh,mlw,w,s,
      b.surface.getSpreadRate(SpeedUnits::FeetPerMinute),
      b.surface.getFirelineIntensity(FirelineIntensityUnits::BtusPerFootPerSecond),
      b.surface.getFlameLength(LengthUnits::Feet),
      b.surface.getReactionIntensity(HeatSourceAndReactionIntensityUnits::BtusPerSquareFootPerMinute));
  }}
}
