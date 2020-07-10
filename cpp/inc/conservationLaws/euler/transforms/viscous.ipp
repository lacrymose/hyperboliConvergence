
// set <-> state

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
      qv[nd]  =state.temperature();
      qv[nd+1]=state.pressure();

      return qv;
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

// delta <-> delta

   template<EulerViscousDelta   ViscDelT,
            EulerConservedDelta ConsDelT,
            EulerState            StateT,
            floating_point          Real>
      requires ConsistentTypes<LawType::Euler,dim_of_v<StateT>,Real,
                               ViscDelT,ConsDelT,StateT>
   ViscDelT delta2Delta( const Species<LawType::Euler,Real>& species,
                         const StateT&                         state,
                         const ConsDelT&                         dqc )
  {
      constexpr int nDim = dim_of_v<StateT>;

      ViscDelT dqv;

      const Real dr  = dqc[nDim];
      const Real dre = dqc[nDim+1];

      const Real rho1 = 1./state.density();
      const Real k = 0.5*state.velocity2();

   // intermediate value used for both temperature and pressure
      Real c = k*dr + dre;

   // velocities
      for( int i=0; i<nDim; i++ )
     {
         const Real   u = state.velocity(i);
         const Real dru = dqc[i];
         const Real  du = ( dru - u*dr )*rho1;

         dqv[i] = du;

         c-=u*dru;
     }

   // pressure
      const Real dp = c*( species.gamma - 1. );

   // temperature
      const Real cv = species.R*species.gamma1;

      const Real t = state.temperature();

      const Real dt = ( c - cv*t*dr )/( cv*state.density() );

      dqv[nDim  ]=dt;
      dqv[nDim+1]=dp;

      return dqv;
  }
   
   template<EulerConservedDelta ConsDelT,
            EulerViscousDelta   ViscDelT,
            EulerState            StateT,
            floating_point          Real>
      requires ConsistentTypes<LawType::Euler,dim_of_v<StateT>,Real,
                               ConsDelT,ViscDelT,StateT>
   ConsDelT delta2Delta( const Species<LawType::Euler,Real>& species,
                         const StateT&                         state,
                         const ViscDelT&                         dqv )
  {
      constexpr int nDim = dim_of_v<StateT>;

      ConsDelT dqc;

      const Real dt = dqv[nDim];
      const Real dp = dqv[nDim+1];

      const Real r = state.density();
      const Real t = state.temperature();
      const Real p = state.pressure();
      const Real k = 0.5*state.velocity2();

   // density
      const Real dr = r*( dp/p - dt/t );

   // drk = k*dr + r*dk needed for dre
      Real drk = k*dr;

   // momentum
      for( int i=0; i<nDim; i++ )
     {
         const Real   u = state.velocity(i);
         const Real  du = dqv[i];
         const Real dru = r*du + u*dr;

         dqc[i] = dru;

         drk+= r*u*du;
     }

   //  re =  p/(gamma-1) +  rk
   // dre = dp/(gamma-1) + drk
      const Real dre = drk + dp*species.gamma1;

      dqc[nDim  ] = dr;
      dqc[nDim+1] = dre;

      return dqc;
  }
   

