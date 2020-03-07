
namespace IdealGas2D
{
   template<>
   template<>
   inline VariableDelta<Conserved>::VariableDelta( const Species& gas, const State& state, const VariableDelta<Conserved>& dq0 )
  {
      var[0]=dq0[0];
      var[1]=dq0[1];
      var[2]=dq0[2];
      var[3]=dq0[3];
  }

   template<>
   template<>
   inline VariableDelta<Conserved>::VariableDelta( const Species& gas, const State& state, const VariableDelta<Viscous>& dq0 )
  {
      Types::Real u,v;
      Types::Real r,h,t,a2;
      Types::Real cp;

      Types::Real dr,dru,drv,dre;
      Types::Real du,dv,dt,dp;

      Types::Real a1,      a4;
      Types::Real b1,b2,   b4;
      Types::Real c1,   c3,c4;
      Types::Real d1,d2,d3,d4;

      du = dq0[0];
      dv = dq0[1];
      dt = dq0[2];
      dp = dq0[3];

      u = state.velocityX();
      v = state.velocityY();
      r = state.density();
      h = state.specificTotalEnthalpy();
      t = state.temperature();
      a2= state.speedOfSound2();

      cp = gas.gamma*gas.Rgas;
      cp/= gas.gamma - 1;

   // construct jacobian
      a1 = gas.gamma/a2;
      b1 = a1*u;
      c1 = a1*v;
      d1 = a1*h - 1;

      b2 = r;
      c3 = r;

      d2 = r*u;
      d3 = r*v;

      a4 = -r/t;
      b4 = a4*u;
      c4 = a4*v;
      d4 = a4*h + r*cp;

   // Ax=b

      dr = a1*dp +                 a4*dt;
      dru= b1*dp + b2*du +         b4*dt;
      drv= c1*dp +         c3*dv + c4*dt;
      dre= d1*dp + d2*du + d3*dv + d4*dt;

      var[0] = dr;
      var[1] = dru;
      var[2] = drv;
      var[3] = dre;

      return;
  }
}
