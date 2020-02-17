
namespace IdealGas2D
{
   inline State::State( Species& gas, ConservedVariables& qc )
  {
      float p,r,t,h,v2,a2;
      float u,v;
      float ru,rv,re;

      float r1,k,gam1;

      gam1=gas.gamma-1;
      gam1=1./gam1;

      r = qc[0];
      ru= qc[1];
      rv= qc[2];
      re= qc[3];

      r1 =  1./r;
      u  = ru*r1;
      v  = rv*r1;

      v2 = u*u;
      v2+= v*v;

      k = 0.5*v2;

      p = ( gas.gamma-1 )*( re-r*k );

      t = p / ( gas.Rgas*r );

      a2 = gas.gamma*gas.Rgas*t;

      h = a2*gam1;
      h+= k;

      state[0] = u;
      state[1] = v;
      state[2] = p;
      state[3] = r;
      state[4] = t;
      state[5] = h;
      state[6] = v2;
      state[7] = a2;
  }

   inline State::State( Species& gas, ViscousVariables& qv )
  {
      float p,r,t,h,v2,a2;
      float u,v;
      float k,gam1;

      gam1=gas.gamma-1;
      gam1=1./gam1;

      u = qv[0];
      v = qv[1];
      t = qv[2];
      p = qv[3];

      r = p/( gas.Rgas*t );

      v2 = u*u;
      v2+= v*v;

      k = 0.5*v2;

      a2 = gas.gamma*gas.Rgas*t;

      h = a2*gam1;
      h+= k;

      state[0] = u;
      state[1] = v;
      state[2] = p;
      state[3] = r;
      state[4] = t;
      state[5] = h;
      state[6] = v2;
      state[7] = a2;
  }
}
