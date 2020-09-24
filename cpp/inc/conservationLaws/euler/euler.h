
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>

// ---------- integral values ----------

   template<> constexpr int nScalarQuantities< LawType::Euler > = 2;
   template<> constexpr int nVectorQuantities< LawType::Euler > = 1;

/*
 * available bases for Variable Sets and Deltas are stored as:
 *    Conserved: velocity, scalar
 */
   enum struct EulerBases
  {
      Conserved,             // { momentum, density,      total energy }
      Characteristic,        // { entropy,  l/r acoustic, vorticity    }
   // Entropy,               // { velocity, entropy,      pressure     }
      Primitive,             // { velocity, density,      pressure     }
      Viscous                // { velocity, temperature,  pressure     }
  };

   template<>
   struct BasisTypeHelper<LawType::Euler> { using type = EulerBases; };

/*
 * available boundary conditions for EulerEquations
 */
   enum struct EulerBCs
  {
      Periodic,
      Riemann,
      InviscidWall
  };

   template<>
   struct BoundaryTypeHelper<LawType::Euler> { using type = EulerBCs; };

// ---------- Law specific types ----------

   template<floating_point Real>
   struct Species<LawType::Euler,Real>
  {
      Real gamma;   // ratio of specific heats
      Real  minf;   // background mach number
      Real  lref;   // lengthscale of the domain
      Real    nu;   // kinematic viscosity
      Real    pr;   // prandtl number
      Real    dt;   // simulation timestep size
      Real     R;   // gas constant

      Real gamma1;   // 1. / ( gamma-1. )
  };

   template<int nDim, floating_point Real>
   struct State<LawType::Euler,nDim,Real>
  {
      std::array<Real,
                 nDim+6> state{0};

      Real& pressure()              { return state[0]; }
      Real& density()               { return state[1]; }
      Real& temperature()           { return state[2]; }
      Real& specificTotalEnthalpy() { return state[3]; }
      Real& velocity2()             { return state[4]; }
      Real& speedOfSound2()         { return state[5]; }
      Real& velocity( const int i ) { return state[6+i]; }

      const Real& pressure()              const { return state[0]; }
      const Real& density()               const { return state[1]; }
      const Real& temperature()           const { return state[2]; }
      const Real& specificTotalEnthalpy() const { return state[3]; }
      const Real& velocity2()             const { return state[4]; }
      const Real& speedOfSound2()         const { return state[5]; }
      const Real& velocity( const int i ) const { return state[6+i]; }
  };


// ---------- Law specific Concepts ----------

   template<typename T>
   concept bool EulerState = 
      is_State_v<T>
   && has_same_law_v<T,law_constant<LawType::Euler>>;

   template<typename T>
   concept bool EulerVarSet =
      has_same_law_v<T,law_constant<LawType::Euler>>
   && ImplementedVarSet<T>;

   template<typename T>
   concept bool EulerVarDelta =
      has_same_law_v<T,law_constant<LawType::Euler>>
   && ImplementedVarDelta<T>;

   // conserved variables

   template<typename T>
   concept bool EulerConservedVariables =
      is_specialised_VarSet_v<T,
                              LawType::Euler,
                              EulerBases::Conserved>;

   template<typename T>
   concept bool EulerConservedDelta =
      is_specialised_VarDelta_v<T,
                                LawType::Euler,
                                EulerBases::Conserved>;

   // primitive variables

   template<typename T>
   concept bool EulerPrimitiveVariables =
      is_specialised_VarSet_v<T,
                              LawType::Euler,
                              EulerBases::Primitive>;

   template<typename T>
   concept bool EulerPrimitiveDelta =
      is_specialised_VarDelta_v<T,
                                LawType::Euler,
                                EulerBases::Primitive>;

   // viscous variables

   template<typename T>
   concept bool EulerViscousVariables =
      is_specialised_VarSet_v<T,
                              LawType::Euler,
                              EulerBases::Viscous>;

   template<typename T>
   concept bool EulerViscousDelta =
      is_specialised_VarDelta_v<T,
                                LawType::Euler,
                                EulerBases::Viscous>;

   // characteristic deltas

   template<typename T>
   concept bool EulerCharacteristicDelta =
      is_specialised_VarDelta_v<T,
                                LawType::Euler,
                                EulerBases::Characteristic>;


// ---------- exact physical flux and spectral radius ----------

   template<int nDim, floating_point Real>
   FluxResult<LawType::Euler,nDim,Real> exactFlux( const Species<LawType::Euler,Real>&  species,
                                                   const geom::Direction<nDim,Real>&     normal,
                                                   const State<LawType::Euler,nDim,Real>& state );

   template<int nDim, floating_point Real>
   Real spectralRadius( const geom::Direction<nDim,Real>&     normal,
                        const State<LawType::Euler,nDim,Real>& state )
  {
      return std::abs( projectedVelocity( normal, state ) ) + std::sqrt( state.speedOfSound2() );
  }

   template<int nDim, floating_point Real>
   WaveSpeeds<LawType::Euler,nDim,Real> wavespeeds( const Species<LawType::Euler,Real>&  species,
                                                    const geom::Direction<nDim,Real>&     normal,
                                                    const State<LawType::Euler,nDim,Real>& state )
  {
      WaveSpeeds<LawType::Euler,nDim,Real> speeds;
      const Real u = projectedVelocity( normal, state );
      const Real a = sqrt(state.speedOfSound2());
      speeds[0]=u;
      speeds[1]=u-a;
      speeds[2]=u+a;
      for( unsigned int i=0; i<nDim-1; i++ ){ speeds[3+i]=u; }
      return speeds;
  }


// ---------- transformation functions set<->state ----------

// conserved variables

   template<EulerConservedVariables ConsVarT, EulerState StateT, floating_point Real>
      requires   SameDim<   ConsVarT,StateT>
              && SameFPType<ConsVarT,StateT,Real>
   ConsVarT state2Set( const Species<LawType::Euler,Real>& species, const StateT& state );

   template<EulerConservedVariables ConsVarT, floating_point Real>
      requires SameFPType<ConsVarT,Real>
   state_t<ConsVarT> set2State( const Species<LawType::Euler,Real>& species, const ConsVarT& qc );


// primitive variables

   template<EulerPrimitiveVariables PrimVarT, EulerState StateT, floating_point Real>
      requires   SameDim<   PrimVarT,StateT>
              && SameFPType<PrimVarT,StateT,Real>
   PrimVarT state2Set( const Species<LawType::Euler,Real>& species, const StateT& state );

   template<EulerPrimitiveVariables PrimVarT, floating_point Real>
      requires SameFPType<PrimVarT,Real>
   state_t<PrimVarT> set2State( const Species<LawType::Euler,Real>& species, const PrimVarT& qp );


// viscous variables

   template<EulerViscousVariables ViscVarT, EulerState StateT, floating_point Real>
      requires   SameDim<   ViscVarT,StateT>
              && SameFPType<ViscVarT,StateT,Real>
   ViscVarT state2Set( const Species<LawType::Euler,Real>& species, const StateT& state );

   template<EulerViscousVariables ViscVarT, floating_point Real>
      requires SameFPType<ViscVarT,Real>
   state_t<ViscVarT> set2State( const Species<LawType::Euler,Real>& species, const ViscVarT& qv );


// ---------- transformation functions delta<->delta ----------

// conserved <-> viscous

   // conserved -> viscous
   template<EulerViscousDelta   ViscDelT,
            EulerConservedDelta ConsDelT,
            EulerState            StateT,
            floating_point          Real>
      requires ConsistentTypes<LawType::Euler,dim_of_v<StateT>,Real,
                               ViscDelT,ConsDelT,StateT>
   ViscDelT delta2Delta( const Species<LawType::Euler,Real>& species,
                         const StateT&                         state,
                         const ConsDelT&                         dqc );
   
   // viscous -> conserved
   template<EulerConservedDelta ConsDelT,
            EulerViscousDelta   ViscDelT,
            EulerState            StateT,
            floating_point          Real>
      requires ConsistentTypes<LawType::Euler,dim_of_v<StateT>,Real,
                               ConsDelT,ViscDelT,StateT>
   ConsDelT delta2Delta( const Species<LawType::Euler,Real>& species,
                         const StateT&                         state,
                         const ViscDelT&                         dqv );
   
// conserved <-> primitive

   // conserved -> primitive
   template<EulerPrimitiveDelta PrimDelT,
            EulerConservedDelta ConsDelT,
            EulerState            StateT,
            floating_point          Real>
      requires ConsistentTypes<LawType::Euler,dim_of_v<StateT>,Real,
                               PrimDelT,ConsDelT,StateT>
   PrimDelT delta2Delta( const Species<LawType::Euler,Real>& species,
                         const StateT&                         state,
                         const ConsDelT&                         dqc );
   
   // primitive -> conserved
   template<EulerConservedDelta ConsDelT,
            EulerPrimitiveDelta PrimDelT,
            EulerState            StateT,
            floating_point          Real>
      requires ConsistentTypes<LawType::Euler,dim_of_v<StateT>,Real,
                               ConsDelT,PrimDelT,StateT>
   ConsDelT delta2Delta( const Species<LawType::Euler,Real>& species,
                         const StateT&                         state,
                         const PrimDelT&                         dqp );

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
                         const CharDelT&                         dqw );
   
   // primitive -> characteristic
   template<EulerCharacteristicDelta CharDelT,
            EulerPrimitiveDelta      PrimDelT,
            EulerState                 StateT,
            floating_point               Real>
      requires ConsistentTypes<LawType::Euler,dim_of_v<StateT>,Real,
                               CharDelT,PrimDelT,StateT>
   CharDelT delta2Delta( const Species<LawType::Euler,Real>& species,
                         const StateT&                         state,
                         const PrimDelT&                         dqp );


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
                         const CharDelT&                         dqw );
   
   // conserved -> characteristic
   template<EulerCharacteristicDelta CharDelT,
            EulerConservedDelta      ConsDelT,
            EulerState                 StateT,
            floating_point               Real>
      requires ConsistentTypes<LawType::Euler,dim_of_v<StateT>,Real,
                               CharDelT,ConsDelT,StateT>
   CharDelT delta2Delta( const Species<LawType::Euler,Real>& species,
                         const StateT&                         state,
                         const ConsDelT&                         dqc );



// ---------- standard ideal gas species ----------

   template<floating_point Real>
   Species<LawType::Euler,Real> get_air_species();


// ---------- flux functions ----------

/*
 * identifier tags for diffusion at low mach number asymptotic solutions
 */
   enum struct LowMachScaling
  {
      Convective,             // scale diffusion to converge to the convective solution with no acoustics (second order pressure)
      Acoustic,               // scale diffusion to converge to the acoustic solution (first order pressure)
      Adaptive                // scale diffusion adaptively based on ratio of timestep and acoustic timescale
  };

/*
 * AUSM+up flux for all speeds - M. S. Liou 2006
 *    includes template options for convective/acoustic/adaptive scaling of pressure/velocity at low Mach number
 */
   template<LowMachScaling VFluxScaling = LowMachScaling::Convective,
            LowMachScaling PFluxScaling = LowMachScaling::Convective>
   struct AusmPlusUP : FluxInterface<AusmPlusUP<VFluxScaling,
                                                PFluxScaling>,
                                     LawType::Euler>
  {
      template<EulerState StateT, int nDim, floating_point Real>
         requires   SameDim<   StateT,dim_constant<nDim>>
                 && SameFPType<StateT,Real>
      FluxResult<LawType::Euler,nDim,Real> flux( const Species<LawType::Euler,Real>& species,
                                                 const geom::Surface<nDim,Real>&        face,
                                                 const StateT&                            sl,
                                                 const StateT&                            sr ) const;
  };

/*
 * SLAU simple low-dissipation ausm flux - Shima & Kitamura 2011
 *    includes template options for convective/acoustic/adaptive scaling of pressure/velocity at low Mach number
 */
   template<LowMachScaling VFluxScaling = LowMachScaling::Convective,
            LowMachScaling PFluxScaling = LowMachScaling::Acoustic>
   struct Slau : FluxInterface<Slau<VFluxScaling,
                                    PFluxScaling>,
                               LawType::Euler>
  {
      template<EulerState StateT, int nDim, floating_point Real>
         requires   SameDim<   StateT,dim_constant<nDim>>
                 && SameFPType<StateT,Real>
      FluxResult<LawType::Euler,nDim,Real> flux( const Species<LawType::Euler,Real>& species,
                                                 const geom::Surface<nDim,Real>&        face,
                                                 const StateT&                            sl,
                                                 const StateT&                            sr ) const;
  };

/*
 * Roe flux utilities
 */
   template<int nDim, floating_point Real>
   State<LawType::Euler,nDim,Real> roeAverage( const Species<LawType::Euler,Real>& species,
                                               const State<LawType::Euler,nDim,Real>&   sl,
                                               const State<LawType::Euler,nDim,Real>&   sr );

   template<int nDim, floating_point Real>
   WaveSpeeds<LawType::Euler,nDim,Real> entropyfix( const WaveSpeeds<LawType::Euler,nDim,Real>& la,
                                                    const WaveSpeeds<LawType::Euler,nDim,Real>& ll,
                                                    const WaveSpeeds<LawType::Euler,nDim,Real>& lr );

// ---------- implementation files  ----------

# include <conservationLaws/euler/fluxes/exactFlux.ipp>
# include <conservationLaws/euler/fluxes/ausmPlusUP.ipp>
# include <conservationLaws/euler/fluxes/slau.ipp>
# include <conservationLaws/euler/fluxes/roe.ipp>

# include <conservationLaws/euler/transforms/conserved.ipp>
# include <conservationLaws/euler/transforms/characteristic.ipp>
# include <conservationLaws/euler/transforms/primitive.ipp>
# include <conservationLaws/euler/transforms/viscous.ipp>

# include <conservationLaws/euler/species.ipp>

