
namespace IdealGas2D
{
   template<>
   template<>
   inline VariableDelta<Conserved>::VariableDelta( const Species& gas, const State& state, const VariableDelta<Viscous>& dq0 )
  {
      float u,v;
      float r,h,t,a2;
      float cp;

      float dr,dru,drv,dre;
      float du,dv,dt,dp;

      float a1,      a4;
      float b1,b2,   b4;
      float c1,   c3,c4;
      float d1,d2,d3,d4;

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
