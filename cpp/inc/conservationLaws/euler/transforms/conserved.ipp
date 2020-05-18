
   template<EulerConservedVariables ConsVarT, EulerState StateT>
      requires SameDim<ConsVarT,StateT>
   ConsVarT state2Set( const Species<LawType::Euler>& species, const StateT& state )
  {
      constexpr int nd = dim_of_v<StateT>;

   // unpack state
      const Types::Real r = state.density();
      const Types::Real p = state.pressure();
      const Types::Real h = state.specificTotalEnthalpy();

   // total energy
      const Types::Real re = r*h - p;

   // assemble conserved variables
      // { momentum, density, total energy }
      ConsVarT qc;

      for( int i=0; i<nd; i++ )
     {
         qc[i] = r*state.velocity(i);
     }
      qc[nd]  =r;
      qc[nd+1]=re;

      return qc;
  }

   template<EulerConservedVariables ConsVarT>
   state_t<ConsVarT> set2State( const Species<LawType::Euler>& species, const ConsVarT& qc )
  {
      using StateT = state_t<ConsVarT>;
      constexpr int nd = dim_of_v<StateT>;

      const Types::Real gam  = species.gamma;
      const Types::Real R    = species.R;
      const Types::Real gam1 = 1./( gam - 1. );

   // unpack conserved variables
      // density and total energy
      const Types::Real r = qc[nd];
      const Types::Real re= qc[nd+1];
      const Types::Real r1= 1./r;

   // velocities
      StateT state;
      state.velocity2()=0;
      for( int i=0; i<nd; i++ )
     {
         state.velocity(i) = qc[i]*r1;
     }
      for( int i=0; i<nd; i++ )
     {
         state.velocity2() += state.velocity(i)*state.velocity(i);
     }
      const Types::Real k = 0.5*state.velocity2();

   // thermodynamic quantities
      const Types::Real p = ( gam-1. )*( re - r*k );
      const Types::Real t = p / ( R*r );
      const Types::Real a2= gam*R*t;
      const Types::Real h = a2*gam1 + k;

   // assemble state

      state.pressure()      = p;
      state.density()       = r;
      state.temperature()   = t;
      state.specificTotalEnthalpy() = h;
      state.speedOfSound2() = a2;

      return state;
  }

