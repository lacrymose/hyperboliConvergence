
namespace IdealGas2D
{
   template<>
   template<>
   inline VariableDelta<Viscous>::VariableDelta( const Species& gas, const State& state, const VariableDelta<Viscous>& dq0 )
  {
      var[0]=dq0[0];
      var[1]=dq0[1];
      var[2]=dq0[2];
      var[3]=dq0[3];
  }

   template<>
   template<>
   inline VariableDelta<Viscous>::VariableDelta( const Species& gas, const State& state, const VariableDelta<Conserved>& dq0 )
  {
      float u,v;
      float r,h,t,a2;
      float cp;

      float dr,dru,drv,dre;
      float du,dv,dt,dp;

      float a1_1;

      float a1,          a4;
      float b1,b2_1,     b4;
      float c1,     c3_1,c4;
      float d1,d2,  d3,  d4;

      float omega1,omega2,omega3;

      dr = dq0[0];
      dru= dq0[1];
      drv= dq0[2];
      dre= dq0[3];

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
      a1_1 = 1./a1;

      b2_1 = 1./r;
      c3_1 = b2_1;

      d2 = r*u;
      d3 = r*v;

      a4 = -r/t;
      b4 = a4*u;
      c4 = a4*v;
      d4 = a4*h + r*cp;

   // solve linear system
      omega1 = d1 - d2*b1*b2_1 - d3*c1*c3_1;
      omega2 = d4 - d2*b4*b2_1 - d3*c4*c3_1;

      omega3 = omega1*dr*a1_1 + d2*dru*b2_1 + d3*drv*c3_1;

      dt = ( dre - omega3 )/ ( omega2 - omega1*a4*a1_1 );

      dp = ( dr - a4*dt )*a1_1;

      du = ( dru - b1*dp - b4*dt )*b2_1;
      dv = ( drv - c1*dp - c4*dt )*c3_1;

      var[0] = du;
      var[1] = dv;
      var[2] = dt;
      var[3] = dp;

      return;
  }
}
