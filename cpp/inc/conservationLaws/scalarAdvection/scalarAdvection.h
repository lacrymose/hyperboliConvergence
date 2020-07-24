
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>

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

   template<floating_point Real>
   struct Species<LawType::ScalarAdvection,Real> {};

   template<int nDim, floating_point Real>
   struct State<LawType::ScalarAdvection,nDim,Real>
  {
      std::array<Real,
                 nVar<LawType::ScalarAdvection,nDim>> state{0};

            Real& scalar()       { return state[0]; }
      const Real& scalar() const { return state[0]; }

            Real& velocity( const int i )       { return state[1+i]; }
      const Real& velocity( const int i ) const { return state[1+i]; }
  };


// ---------- Law specific Concepts ----------

   template<typename T>
   concept bool ScalarState = 
      is_State_v<T>
   && has_same_law_v<T,law_constant<LawType::ScalarAdvection>>;

   template<typename T>
   concept bool ScalarVarSet =
      has_same_law_v<T,law_constant<LawType::ScalarAdvection>>
   && ImplementedVarSet<T>;

   template<typename T>
   concept bool ScalarConservedVariables =
      is_specialised_VarSet_v<T,
                              LawType::ScalarAdvection,
                              ScalarAdvectionBases::Conserved>;


// ---------- exact physical flux and spectral radius ----------

   template<int nDim, floating_point Real>
   FluxResult<LawType::ScalarAdvection,nDim,Real> exactFlux( const Species<LawType::ScalarAdvection,Real>&  species,
                                                             const geom::Direction<nDim,Real>&               normal,
                                                             const State<LawType::ScalarAdvection,nDim,Real>& state );

   template<int nDim, floating_point Real>
   Real spectralRadius( const geom::Direction<nDim,Real>&               normal,
                        const State<LawType::ScalarAdvection,nDim,Real>& state )
  {
      return std::abs( projectedVelocity( normal, state ) );
  }


// ---------- transformation functions ----------

   template<ScalarConservedVariables ConsVarT, ScalarState StateT, floating_point Real>
      requires   SameDim<ConsVarT,StateT>
              && SameFPType<ConsVarT,StateT,Real>
   ConsVarT state2Set( const Species<LawType::ScalarAdvection,Real>& species, const StateT& state );

   template<ScalarConservedVariables ConsVarT, floating_point Real>
      requires SameFPType<ConsVarT,Real>
   state_t<ConsVarT> set2State( const Species<LawType::ScalarAdvection,Real>& species, const ConsVarT& qc );


# include <conservationLaws/scalarAdvection/fluxes/exactFlux.ipp>

# include <conservationLaws/scalarAdvection/transforms/conserved.ipp>
