
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

   template< char C >
   inline State::State( const Species& gas, const VariableSet<C>& q )
  {
      std::cout << std::endl;
      std::cout << "Warning:" << std::endl;
      std::cout << "State::State( const Species& gas, const VariableSet<D>& q0 )"
                << " is not yet defined for "
                << "C = '" << C << "'" << std::endl;
      std::cout << std::endl;
      assert( false );
  }

   template<>
   inline State::State( const Species& gas, const VariableSet<'c'>& q )
  {
      float p,r,t,h,v2,a2;
      float u,v;
      float ru,rv,re;

      float r1,k,gam1;

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
   inline State::State( const Species& gas, const VariableSet<'v'>& q )
  {
      float p,r,t,h,v2,a2;
      float u,v;
      float k,gam1;

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
}
