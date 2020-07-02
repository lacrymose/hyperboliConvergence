
   template<EulerPrimitiveVariables PrimVarT, EulerState StateT, floating_point Real>
      requires   SameDim<PrimVarT,StateT>
              && SameFPType<PrimVarT,StateT,Real>
   PrimVarT state2Set( const Species<LawType::Euler,Real>& species, const StateT& state )
  {
      PrimVarT qp;

      constexpr int nd = dim_of_v<StateT>;
      for( int i=0; i<nd; i++ )
     {
         qp[i] = state.velocity(i);
     }
      qp[nd]  =state.density();
      qp[nd+1]=state.pressure();

      return qp;
  }

   template<EulerPrimitiveVariables PrimVarT, floating_point Real>
      requires SameFPType<PrimVarT,Real>
   state_t<PrimVarT> set2State( const Species<LawType::Euler,Real>& species, const PrimVarT& qp )
  {
      using StateT = state_t<PrimVarT>;
      constexpr int nd = dim_of_v<PrimVarT>;

      const Real gam  = species.gamma;
      const Real R    = species.R;
      const Real gam1 = species.gamma1;

      StateT state;

   // velocities
      state.velocity2()=0;
      for( int i=0; i<nd; i++ )
     {
         state.velocity(i) = qp[i];
         state.velocity2()+= qp[i]*qp[i];
     }
      const Real k = 0.5*state.velocity2();

   // thermodynamic quantities
      const Real r = qp[nd];
      const Real p = qp[nd+1];
      const Real t = p/( r*R );
      const Real a2= gam*R*t;
      const Real h = a2*gam1 + k;

      state.specificTotalEnthalpy() = h;
      state.pressure()      = p;
      state.density()       = r;
      state.temperature()   = t;
      state.speedOfSound2() = a2;

      return state;
  }

