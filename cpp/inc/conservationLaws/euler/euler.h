
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <types.h>

// ---------- integral values ----------

/*
 * available bases for Variable Sets and Deltas are stored as:
 *    Conserved: velocity, scalar
 */
   enum struct EulerBases { Conserved,       // { momentum, density,     total energy }
                            Characteristic,  // { entropy,  vorticity,   left/right acoustic }
                            Viscous,         // { velocity, temperature, pressure }
                            Primitive };     // { velocity, density,     pressure }

   template<>
   struct BasisTypeHelper<LawType::Euler> { using type = EulerBases; };

   template<> constexpr int nScalarQuantities< LawType::Euler > = 2;
   template<> constexpr int nVectorQuantities< LawType::Euler > = 1;


// ---------- Law specific types ----------

   template<>
   struct Species<LawType::Euler>
  {
      Types::Real gamma;   // ratio of specific heats
      Types::Real  minf;   // background mach number
      Types::Real    nu;   // kinematic viscosity
      Types::Real    pr;   // prandtl number
      Types::Real    dt;   // simulation timescale
      Types::Real     R;   // gas constant
  };

   template<int nDim>
   struct State<LawType::Euler,nDim>
  {
      std::array<Types::Real,
                 nDim+6> state{0};

      Types::Real& pressure()              { return state[0]; }
      Types::Real& density()               { return state[1]; }
      Types::Real& temperature()           { return state[2]; }
      Types::Real& specificTotalEnthalpy() { return state[3]; }
      Types::Real& velocity2()             { return state[4]; }
      Types::Real& speedOfSound2()         { return state[5]; }
      Types::Real& velocity( const int i ) { return state[6+i]; }

      const Types::Real& pressure()              const { return state[0]; }
      const Types::Real& density()               const { return state[1]; }
      const Types::Real& temperature()           const { return state[2]; }
      const Types::Real& specificTotalEnthalpy() const { return state[3]; }
      const Types::Real& velocity2()             const { return state[4]; }
      const Types::Real& speedOfSound2()         const { return state[5]; }
      const Types::Real& velocity( const int i ) const { return state[6+i]; }
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

   template<int nDim>
   FluxResult<LawType::Euler,nDim> exactFlux( const Species<LawType::Euler>&   species,
                                              const Geometry::Direction<nDim>&  normal,
                                              const State<LawType::Euler,nDim>&  state );


// ---------- transformation functions ----------

// conserved variables
   template<EulerConservedVariables ConsVarT, EulerState StateT>
      requires SameDim<ConsVarT,StateT>
   ConsVarT state2Set( const Species<LawType::Euler>& species, const StateT& state );

   template<EulerConservedVariables ConsVarT>
   state_t<ConsVarT> set2State( const Species<LawType::Euler>& species, const ConsVarT& qc );

// viscous variables
   template<EulerViscousVariables ViscVarT, EulerState StateT>
      requires SameDim<ViscVarT,StateT>
   ViscVarT state2Set( const Species<LawType::Euler>& species, const StateT& state );

   template<EulerViscousVariables ViscVarT>
   state_t<ViscVarT> set2State( const Species<LawType::Euler>& species, const ViscVarT& qv );



# include <conservationLaws/euler/fluxes/exactFlux.ipp>

# include <conservationLaws/euler/transforms/conserved.ipp>
# include <conservationLaws/euler/transforms/viscous.ipp>

# include <conservationLaws/euler/species.ipp>

