
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <types.h>

// ---------- integral values ----------

/*
 * available bases for Variable Sets and Deltas are stored as:
 *    Conserved: velocity, scalar
 */
   enum struct ScalarAdvectionBases { Conserved,
                                      Characteristic };

   template<>
   struct BasisTypeHelper<LawType::ScalarAdvection> { using type = ScalarAdvectionBases; };

   template<> constexpr int nScalarQuantities< LawType::ScalarAdvection > = 1;
   template<> constexpr int nVectorQuantities< LawType::ScalarAdvection > = 1;


// ---------- Law specific types ----------

   template<>
   struct Species<LawType::ScalarAdvection>
  {
      Types::Real scale=1;
  };

   template<int nDim>
   struct State<LawType::ScalarAdvection,nDim>
  {
      std::array<Types::Real,
                 nVar<LawType::ScalarAdvection,nDim>> state{0};

            Types::Real& scalar()       { return state[0]; }
      const Types::Real& scalar() const { return state[0]; }

            Types::Real& velocity( const int i )       { return state[1+i]; }
      const Types::Real& velocity( const int i ) const { return state[1+i]; }
  };


// ---------- Law specific Concepts ----------

   template<typename T>
   concept bool ScalarState = 
      is_State_v<T>
   && has_same_law_v<T,law_constant<LawType::ScalarAdvection>>;

   template<typename T>
   concept bool ScalarConservedVariables =
//    ImplementedVarSet<T>
      is_specialised_VarSet_v<T,
                              LawType::ScalarAdvection,
                              ScalarAdvectionBases::Conserved>;


// ---------- exact physical flux ----------

   template<int nDim>
   FluxResult<LawType::ScalarAdvection,nDim> exactFlux( const Species<LawType::ScalarAdvection>&  species,
                                                        const Geometry::Direction<nDim>&           normal,
                                                        const State<LawType::ScalarAdvection,nDim>& state );


// ---------- transformation functions ----------

   template<ScalarConservedVariables ConsVarT, ScalarState StateT>
      requires SameDim<ConsVarT,StateT>
   ConsVarT state2Set( const Species<LawType::ScalarAdvection>& species, const StateT& state );

   template<ScalarConservedVariables ConsVarT>
   state_t<ConsVarT> set2State( const Species<LawType::ScalarAdvection>& species, const ConsVarT& qc );


# include <conservationLaws/scalarAdvection/fluxes/exactFlux.ipp>

# include <conservationLaws/scalarAdvection/transforms/conserved.ipp>
