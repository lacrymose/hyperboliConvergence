
namespace IdealGas2D
{
   template<>
   template<>
   inline VariableSet<Conserved>::VariableSet( const Species& gas, const VariableSet<Viscous>& q0 )
  {
      float u,v,t,p;
      float r,ru,rv,re;
      float k;
      float gam1;

      gam1=gas.gamma-1;
      gam1=1./gam1;

      u = q0[0];
      v = q0[1];
      t = q0[2];
      p = q0[3];

      r = gas.Rgas*t;
      r = p/r;

      ru = r*u;
      rv = r*v;

      k =   u*u;
      k+=   v*v;
      k*= 0.5*r;

      re = p*gam1;
      re+=      k;

      var[0]=r;
      var[1]=ru;
      var[2]=rv;
      var[3]=re;

      return;
  }

   template<>
   inline VariableSet<Conserved>::VariableSet( const Species& gas, const State& state )
  {
      float r,ru,rv,re;
      float u,v,p,h;

      r = state.density();
      u = state.velocityX();
      v = state.velocityY();
      p = state.pressure();
      h = state.specificTotalEnthalpy();

      ru = r*u;
      rv = r*v;
      re = r*h - p;

      var[0]=r;
      var[1]=ru;
      var[2]=rv;
      var[3]=re;

      return;
  }
}
