
// primitive <-> characteristic

   // characteristic -> primitive
   template<EulerPrimitiveDelta      PrimDelT,
            EulerCharacteristicDelta CharDelT,
            EulerState                 StateT,
            floating_point               Real>
      requires ConsistentTypes<LawType::Euler,dim_of_v<StateT>,Real,
                               PrimDelT,CharDelT,StateT>
   PrimDelT delta2Delta( const Species<LawType::Euler,Real>& species,
                         const StateT&                         state,
                         const CharDelT&                         dqw )
  {
      constexpr int nDim = dim_of_v<StateT>;

      const Real r = state.density();
      const Real a = sqrt( state.speedOfSound2() );

      PrimDelT dqp;

   // pressure
      const Real dp = 0.5*r*a*( dqw[2] - dqw[1] );

   // density
      const Real dr = dp/(a*a) - r*dqw[0]/species.gamma;

   // normal velocity
      dqp[0] = 0.5*( dqw[1] + dqw[2] );

   // vorticity
      for( int i=1; i<nDim; i++ )
     {
         dqp[i] = dqw[2+i];
     }

      dqp[nDim  ] = dr;
      dqp[nDim+1] = dp;

      return dqp;
  }
   
   // primitive -> characteristic
   template<EulerCharacteristicDelta CharDelT,
            EulerPrimitiveDelta      PrimDelT,
            EulerState                 StateT,
            floating_point               Real>
      requires ConsistentTypes<LawType::Euler,dim_of_v<StateT>,Real,
                               CharDelT,PrimDelT,StateT>
   CharDelT delta2Delta( const Species<LawType::Euler,Real>& species,
                         const StateT&                         state,
                         const PrimDelT&                         dqp )
  {
      constexpr int nDim = dim_of_v<StateT>;

      const Real r1= 1./state.density();
      const Real a1= 1./sqrt( state.speedOfSound2() );

      const Real dr = dqp[nDim  ];
      const Real dp = dqp[nDim+1];
      const Real du = dqp[0];

      CharDelT dqw;

   // entropy
      dqw[0] = species.gamma*( dp*a1*a1 - dr )*r1;

   // acoustics
      dqw[1] = du - dp*r1*a1;
      dqw[2] = du + dp*r1*a1;

   // vorticity
      for( int i=1; i<nDim; i++ )
     {
         dqw[2+i] = dqp[i];
     }

      return dqw;
  }


// conserved <-> characteristic

   // characteristic -> conserved
   template<EulerConservedDelta      ConsDelT,
            EulerCharacteristicDelta CharDelT,
            EulerState                 StateT,
            floating_point               Real>
      requires ConsistentTypes<LawType::Euler,dim_of_v<StateT>,Real,
                               ConsDelT,CharDelT,StateT>
   ConsDelT delta2Delta( const Species<LawType::Euler,Real>& species,
                         const StateT&                         state,
                         const CharDelT&                         dqw )
  {
      constexpr int nDim = dim_of_v<StateT>;
      using PrimDelT = VariableDelta<LawType::Euler,nDim,EulerBases::Primitive,Real>;

      return delta2Delta<ConsDelT>( species,
                                    state,
                                    delta2Delta<PrimDelT>( species,
                                                           state,
                                                           dqw ) );
  }
   
   // conserved -> characteristic
   template<EulerCharacteristicDelta CharDelT,
            EulerConservedDelta      ConsDelT,
            EulerState                 StateT,
            floating_point               Real>
      requires ConsistentTypes<LawType::Euler,dim_of_v<StateT>,Real,
                               CharDelT,ConsDelT,StateT>
   CharDelT delta2Delta( const Species<LawType::Euler,Real>& species,
                         const StateT&                         state,
                         const ConsDelT&                         dqc )
  {
      constexpr int nDim = dim_of_v<StateT>;
      using PrimDelT = VariableDelta<LawType::Euler,nDim,EulerBases::Primitive,Real>;

      return delta2Delta<CharDelT>( species,
                                    state,
                                    delta2Delta<PrimDelT>( species,
                                                           state,
                                                           dqc ) );
  }



