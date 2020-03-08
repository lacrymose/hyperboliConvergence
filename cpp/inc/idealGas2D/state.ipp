
namespace IdealGas2D
{
   inline State::State()
  {
      state[0]=0.;
      state[1]=0.;
      state[2]=0.;
      state[3]=0.;
      state[4]=0.;
      state[5]=0.;
      state[6]=0.;
      state[7]=0.;
  }

   template<>
   inline State::State( const Species& gas, const VariableSet<Conserved>& q )
  {
      Types::Real p,r,t,h,v2,a2;
      Types::Real u,v;
      Types::Real ru,rv,re;

      Types::Real r1,k,gam1;

      gam1=gas.gamma-1;
      gam1=1./gam1;

      r = q[0];
      ru= q[1];
      rv= q[2];
      re= q[3];

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

      return;
  }

   template<>
   inline State::State( const Species& gas, const VariableSet<Viscous>& q )
  {
      Types::Real p,r,t,h,v2,a2;
      Types::Real u,v;
      Types::Real k,gam1;

      gam1=gas.gamma-1;
      gam1=1./gam1;

      u = q[0];
      v = q[1];
      t = q[2];
      p = q[3];

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

      return;
  }

   inline bool operator==( const State& s0, const State& s1 )
  {
      bool same=true;

      same = same and ( fabs( s0.velocityX()             - s1.velocityX()             ) < Types::EPS );
      same = same and ( fabs( s0.velocityY()             - s1.velocityY()             ) < Types::EPS );
      same = same and ( fabs( s0.pressure()              - s1.pressure()              ) < Types::EPS );
      same = same and ( fabs( s0.density()               - s1.density()               ) < Types::EPS );
      same = same and ( fabs( s0.temperature()           - s1.temperature()           ) < Types::EPS );
      same = same and ( fabs( s0.specificTotalEnthalpy() - s1.specificTotalEnthalpy() ) < Types::EPS );
      same = same and ( fabs( s0.velocity2()             - s1.velocity2()             ) < Types::EPS );
      same = same and ( fabs( s0.speedOfSound2()         - s1.speedOfSound2()         ) < Types::EPS );

      return same;
  }

   inline bool operator!=( const State& s0, const State& s1 )
  {
      return !(s0==s1);
  }
}
