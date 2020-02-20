
namespace IdealGas2D
{
   template<>
   template<>
   inline VariableSet<Viscous>::VariableSet( const Species& gas, const VariableSet<Conserved>& q0 )
  {
      float u,v,t,p;
      float r,ru,rv,re;
      float k,r1;
      float gam1;

      gam1=gas.gamma-1;
      gam1=1./gam1;

      r = q0[0];
      ru= q0[1];
      rv= q0[2];
      re= q0[3];

      r1 = 1./r;
      u = ru*r1;
      v = rv*r1;

      k =   u*u;
      k+=   v*v;
      k*= 0.5*r;

      re-=    k;
      p = gas.gamma-1;
      p*=    re;

      t = gas.Rgas*r;
      t = p/t;

      var[0]=u;
      var[1]=v;
      var[2]=t;
      var[3]=p;

      return;
  }

   template<>
   inline VariableSet<Viscous>::VariableSet( const Species& gas, const State& state )
  {
      var[0]=state.velocityX();
      var[1]=state.velocityY();
      var[2]=state.temperature();
      var[3]=state.pressure();

      return;
  }
}
