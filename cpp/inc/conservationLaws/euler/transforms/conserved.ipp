
   template<EulerConservedVariables ConsVarT, EulerState StateT, floating_point Real>
      requires   SameDim<ConsVarT,StateT>
              && SameFPType<ConsVarT,StateT,Real>
   ConsVarT state2Set( const Species<LawType::Euler,Real>& species, const StateT& state )
  {
      constexpr int nd = dim_of_v<StateT>;

   // unpack state
      const Real r = state.density();
      const Real p = state.pressure();
      const Real h = state.specificTotalEnthalpy();

   // total energy
      const Real re = r*h - p;

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

   template<EulerConservedVariables ConsVarT, floating_point Real>
      requires SameFPType<ConsVarT,Real>
   state_t<ConsVarT> set2State( const Species<LawType::Euler,Real>& species, const ConsVarT& qc )
  {
      using StateT = state_t<ConsVarT>;
      constexpr int nd = dim_of_v<StateT>;

      const Real gam  = species.gamma;
      const Real R    = species.R;
      const Real gam1 = species.gamma1;

   // unpack conserved variables
      // density and total energy
      const Real r = qc[nd];
      const Real re= qc[nd+1];
      const Real r1= 1./r;

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
      const Real k = 0.5*state.velocity2();

   // thermodynamic quantities
      const Real p = ( gam-1. )*( re - r*k );
      const Real t = p / ( R*r );
      const Real a2= gam*R*t;
      const Real h = a2*gam1 + k;

   // assemble state

      state.pressure()      = p;
      state.density()       = r;
      state.temperature()   = t;
      state.specificTotalEnthalpy() = h;
      state.speedOfSound2() = a2;

      return state;
  }

