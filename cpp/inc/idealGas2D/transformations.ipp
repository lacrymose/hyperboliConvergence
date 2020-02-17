
namespace IdealGas2D
{
   inline ConservedVariables conservedVariables( Species& gas, ConservedVariables& qc )
  {
      ConservedVariables qc1(qc);
      return qc1;
  }

   inline ConservedVariables conservedVariables( Species& gas, ViscousVariables&   qv )
  {
      ConservedVariables qc;

      float u,v,t,p;
      float r,ru,rv,re;
      float k;
      float gam1;

      gam1=gas.gamma-1;
      gam1=1./gam1;

      u = qv[0];
      v = qv[1];
      t = qv[2];
      p = qv[3];

      r = gas.Rgas*t;
      r = p/r;

      ru = r*u;
      rv = r*v;

      k =   u*u;
      k+=   v*v;
      k*= 0.5*r;

      re = p*gam1;
      re+=      k;

      qc[0]=r;
      qc[1]=ru;
      qc[2]=rv;
      qc[3]=re;

      return qc;
  }

   inline ViscousVariables viscousVariables( Species& gas, ConservedVariables& qc )
  {
      ViscousVariables qv;

      float u,v,t,p;
      float r,ru,rv,re;
      float k,r1;
      float gam1;

      gam1=gas.gamma-1;
      gam1=1./gam1;

      r = qc[0];
      ru= qc[1];
      rv= qc[2];
      re= qc[3];

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

      qv[0]=u;
      qv[1]=v;
      qv[2]=t;
      qv[3]=p;

      return qv;
  }

   inline ViscousVariables viscousVariables( Species& gas, ViscousVariables&   qv )
  {
      ViscousVariables qv1(qv);
      return qv1;
  }

   inline ViscousVariables dViscousVariables( Species& gas, State& state, ConservedVariables& dqc )
  {
      ViscousVariables dqv;

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

      float u,v;
      float r,h,t,a2;
      float cp;

      float dr,dru,drv,dre;
      float du,dv,dt,dp;

      float a1,      a4;
      float b1,b2,   b4;
      float c1,   c3,c4;
      float d1,d2,d3,d4;

      du = dqc[0];
      dv = dqc[1];
      dt = dqc[2];
      dp = dqc[3];

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
