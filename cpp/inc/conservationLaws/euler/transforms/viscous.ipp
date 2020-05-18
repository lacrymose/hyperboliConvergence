
   template<EulerViscousVariables ViscVarT, EulerState StateT>
      requires SameDim<ViscVarT,StateT>
   ViscVarT state2Set( const Species<LawType::Euler>& species, const StateT& state )
  {
      ViscVarT qv;

      constexpr int nd = dim_of_v<StateT>;
      for( int i=0; i<nd; i++ )
     {
         qv[i] = state.velocity(i);
     }
      qv[nd]=state.temperature();
      qv[nd+1]=state.pressure();
  }

   template<EulerViscousVariables ViscVarT>
   state_t<ViscVarT> set2State( const Species<LawType::Euler>& species, const ViscVarT& qv )
  {
      constexpr int nd = dim_of_v<ViscVarT>;
      using StateT = State<LawType::Euler,nd>;

      const Types::Real gam  = species.gamma;
      const Types::Real R    = species.R;
      const Types::Real gam1 = 1./ ( gam-1. );

   // velocities
      StateT state;
      state.velocity2()=0;
      for( int i=0; i<nd; i++ )
     {
         state.velocity(i) = qv[i];
     }
      for( int i=0; i<nd; i++ )
     {
         state.velocity2()+= qv[i]*qv[i];
     }
      const Types::Real k = 0.5*state.velocity2();

   // thermodynamic quantities
      const Types::Real t = qv[nd];
      const Types::Real p = qv[nd+1];
      const Types::Real r = p/( species.R*t );
      const Types::Real a2= gam*R*t;
      const Types::Real h = a2*gam1 + k;

      state.pressure()      = p;
      state.density()       = r;
      state.temperature()   = t;
      state.specificTotalEnthalpy() = h;
      state.speedOfSound2() = a2;

      return state;

  }

