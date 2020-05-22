
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>

// ---------- integral values ----------

/*
 * available bases for Variable Sets and Deltas are stored as:
 *    Conserved: velocity, scalar
 */
   enum struct EulerBases { Conserved,       // { momentum, density,     total energy }
//                          Characteristic,  // { entropy,  vorticity,   left/right acoustic }
//                          Primitive,       // { velocity, density,     pressure }
                            Viscous };       // { velocity, temperature, pressure }

   template<>
   struct BasisTypeHelper<LawType::Euler> { using type = EulerBases; };

   template<> constexpr int nScalarQuantities< LawType::Euler > = 2;
   template<> constexpr int nVectorQuantities< LawType::Euler > = 1;


// ---------- Law specific types ----------

   template<floating_point Real>
   struct Species<LawType::Euler,Real>
  {
      Real gamma;   // ratio of specific heats
      Real  minf;   // background mach number
      Real    nu;   // kinematic viscosity
      Real    pr;   // prandtl number
      Real    dt;   // simulation timescale
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
   concept bool EulerConservedVariables =
      is_specialised_VarSet_v<T,
                              LawType::Euler,
                              EulerBases::Conserved>;
   template<typename T>
   concept bool EulerViscousVariables =
      is_specialised_VarSet_v<T,
                              LawType::Euler,
                              EulerBases::Viscous>;


// ---------- exact physical flux ----------

   template<int nDim, floating_point Real>
   FluxResult<LawType::Euler,nDim,Real> exactFlux( const Species<LawType::Euler,Real>&  species,
                                                   const geom::Direction<nDim,Real>&     normal,
                                                   const State<LawType::Euler,nDim,Real>& state );


// ---------- transformation functions ----------

// conserved variables
   template<EulerConservedVariables ConsVarT, EulerState StateT, floating_point Real>
      requires   SameDim<ConsVarT,StateT>
              && SameFPType<ConsVarT,StateT,Real>
   ConsVarT state2Set( const Species<LawType::Euler,Real>& species, const StateT& state );

   template<EulerConservedVariables ConsVarT, floating_point Real>
      requires SameFPType<ConsVarT,Real>
   state_t<ConsVarT> set2State( const Species<LawType::Euler,Real>& species, const ConsVarT& qc );

// viscous variables
   template<EulerViscousVariables ViscVarT, EulerState StateT, floating_point Real>
      requires   SameDim<ViscVarT,StateT>
              && SameFPType<ViscVarT,StateT,Real>
   ViscVarT state2Set( const Species<LawType::Euler,Real>& species, const StateT& state );

   template<EulerViscousVariables ViscVarT, floating_point Real>
      requires SameFPType<ViscVarT,Real>
   state_t<ViscVarT> set2State( const Species<LawType::Euler,Real>& species, const ViscVarT& qv );

// ---------- standard ideal gas species ----------

   template<floating_point Real>
   Species<LawType::Euler,Real> get_air_species();

# include <conservationLaws/euler/fluxes/exactFlux.ipp>

# include <conservationLaws/euler/transforms/conserved.ipp>
# include <conservationLaws/euler/transforms/viscous.ipp>

# include <conservationLaws/euler/species.ipp>

