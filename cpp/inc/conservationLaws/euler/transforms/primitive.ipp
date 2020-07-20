
// set <-> state

   // state -> set
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

   // set -> state
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

// deltas

   // conserved -> primitive
   template<EulerPrimitiveDelta PrimDelT,
            EulerConservedDelta ConsDelT,
            EulerState            StateT,
            floating_point          Real>
      requires ConsistentTypes<LawType::Euler,dim_of_v<StateT>,Real,
                               PrimDelT,ConsDelT,StateT>
   PrimDelT delta2Delta( const Species<LawType::Euler,Real>& species,
                         const StateT&                         state,
                         const ConsDelT&                         dqc )
  {
      constexpr int nDim = dim_of_v<StateT>;

      PrimDelT dqp;

      const Real dr  = dqc[nDim];
      const Real dre = dqc[nDim+1];

      const Real rho1 = 1./state.density();
      const Real k = 0.5*state.velocity2();

   // intermediate value used for pressure
      Real c = k*dr + dre;

   // velocities
      for( int i=0; i<nDim; i++ )
     {
         const Real   u = state.velocity(i);
         const Real dru = dqc[i];
         const Real  du = ( dru - u*dr )*rho1;

         dqp[i] = du;

         c-=u*dru;
     }

   // pressure
      const Real dp = c*( species.gamma - 1. );

      dqp[nDim  ]=dr;
      dqp[nDim+1]=dp;

      return dqp;
  }
   
   // primitive -> conserved
   template<EulerConservedDelta ConsDelT,
            EulerPrimitiveDelta PrimDelT,
            EulerState            StateT,
            floating_point          Real>
      requires ConsistentTypes<LawType::Euler,dim_of_v<StateT>,Real,
                               ConsDelT,PrimDelT,StateT>
   ConsDelT delta2Delta( const Species<LawType::Euler,Real>& species,
                         const StateT&                         state,
                         const PrimDelT&                         dqp )
  {
      constexpr int nDim = dim_of_v<StateT>;

      ConsDelT dqc;

      const Real dr = dqp[nDim];
      const Real dp = dqp[nDim+1];
      const Real r = state.density();

   // dk = 0.5*( 2udu + 2vdv ... )
      Real dk=0;

   // momentum
      for( int i=0; i<nDim; i++ )
     {
         const Real   u = state.velocity(i);
         const Real  du = dqp[i];
         const Real dru = r*du + u*dr;

         dqc[i] = dru;

         dk+= u*du;
     }

   // total energy
      //  re =  p/(gamma-1) +  rk
      // dre = dp/(gamma-1) + drk
      // dre = dp/(gamma-1) + r*dk + k*dr

      const Real k = 0.5*state.velocity2();

      const Real dre = dp*species.gamma1 + r*dk + k*dr;

      dqc[nDim  ] = dr;
      dqc[nDim+1] = dre;

      return dqc;
  }

