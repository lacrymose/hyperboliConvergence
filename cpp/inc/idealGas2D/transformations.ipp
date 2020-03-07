
namespace IdealGas2D
{
   inline ViscousVariables dViscousVariables( Species& gas, State& state, ConservedVariables& dqc )
  {
      ViscousVariables dqv;

      Types::Real u,v;
      Types::Real r,h,t,a2;
      Types::Real cp;

      Types::Real dr,dru,drv,dre;
      Types::Real du,dv,dt,dp;

      Types::Real a1_1;

      Types::Real a1,          a4;
      Types::Real b1,b2_1,     b4;
      Types::Real c1,     c3_1,c4;
      Types::Real d1,d2,  d3,  d4;

      Types::Real omega1,omega2,omega3;

      dr = dqc[0];
      dru= dqc[1];
      drv= dqc[2];
      dre= dqc[3];

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

      dqv[0] = du;
      dqv[1] = dv;
      dqv[2] = dt;
      dqv[3] = dp;

      return dqv;
  }

   inline ConservedVariables dConservedVariables( Species& gas, State& state, ViscousVariables& dqv )
  {
      ConservedVariables dqc;

      Types::Real u,v;
      Types::Real r,h,t,a2;
      Types::Real cp;

      Types::Real dr,dru,drv,dre;
      Types::Real du,dv,dt,dp;

      Types::Real a1,      a4;
      Types::Real b1,b2,   b4;
      Types::Real c1,   c3,c4;
      Types::Real d1,d2,d3,d4;

      du = dqv[0];
      dv = dqv[1];
      dt = dqv[2];
      dp = dqv[3];

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

      dqc[0] = dr;
      dqc[1] = dru;
      dqc[2] = drv;
      dqc[3] = dre;

      return dqc;
  }
}
