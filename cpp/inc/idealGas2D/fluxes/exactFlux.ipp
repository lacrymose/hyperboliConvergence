
# include <math.h>

namespace IdealGas2D
{
   inline void exactFlux( Species& gas, float n[3], State& s, ConservedVariables& f, float& lmax )
  {
      float r,u,v,h,p,a;
      float un,mn;

   // unpack state
      r = s.density();
      u = s.velocityX();
      v = s.velocityY();
      p = s.pressure();
      h = s.specificTotalEnthalpy();
      a = sqrt( s.speedOfSound2() );

   // face normal velocity and mass flux
      un = u*n[0] + v*n[1];
      mn = r*un;

   // max wavespeed
      lmax = fabs( un ) + a;

   // exact flux
      f[0] = mn;
      f[1] = mn*u + p*n[0];
      f[2] = mn*v + p*n[1];
      f[3] = mn*h;

      return;
  }
}
