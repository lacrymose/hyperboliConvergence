
   template<EulerViscousVariables ViscVarT, EulerState StateT, floating_point Real>
      requires   SameDim<ViscVarT,StateT>
              && SameFPType<ViscVarT,StateT,Real>
   ViscVarT state2Set( const Species<LawType::Euler,Real>& species, const StateT& state )
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

   template<EulerViscousVariables ViscVarT, floating_point Real>
      requires SameFPType<ViscVarT,Real>
   state_t<ViscVarT> set2State( const Species<LawType::Euler,Real>& species, const ViscVarT& qv )
  {
      using StateT = state_t<ViscVarT>;
      constexpr int nd = dim_of_v<ViscVarT>;

      const Real gam  = species.gamma;
      const Real R    = species.R;
      const Real gam1 = species.gamma1;

   // velocities
      StateT state;
      state.velocity2()=0;
      for( int i=0; i<nd; i++ )
     {
         state.velocity(i) = qv[i];
         state.velocity2()+= qv[i]*qv[i];
     }
      const Real k = 0.5*state.velocity2();

   // thermodynamic quantities
      const Real t = qv[nd];
      const Real p = qv[nd+1];
      const Real r = p/( R*t );
      const Real a2= gam*R*t;
      const Real h = a2*gam1 + k;

      state.specificTotalEnthalpy() = h;
      state.pressure()      = p;
      state.density()       = r;
      state.temperature()   = t;
      state.speedOfSound2() = a2;

      return state;

  }

