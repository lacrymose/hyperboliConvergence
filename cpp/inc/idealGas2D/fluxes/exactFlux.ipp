
namespace IdealGas2D
{
   inline void exactFlux( const Species& gas, const float n[3], const State& state, ConservedDelta& f, float& lmax )
  {
      float r,u,v,h,p,a;
      float un,mn;

   // unpack state
      r = state.density();
      u = state.velocityX();
      v = state.velocityY();
      p = state.pressure();
      h = state.specificTotalEnthalpy();
      a = sqrt( state.speedOfSound2() );

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
