
# pragma once

# include <conservationLaws/base/concepts.h>
# include <conservationLaws/base/type-traits.h>
# include <conservationLaws/base/declarations.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>
# include <types.h> 

# include <type_traits>
# include <ostream>

// template<size_t N, typename... Ts>
// concept bool PackOfN = is_pack

// ---------- types for variables in conservation law phase space ----------

/*
 * Vector of variables for a hyperbolic conservation law (Law) in a particular basis for the phase space (Basis) with (nDim) spatial dimensions
 * VariableSet is a point in an affine space, with VariableDelta<Law,nDim,Basis> being displacements in this space
 */
   template<LawType Law, int nDim, BasisType<Law> Basis>
   struct VariableSet
  {
   // equivalent VariableDelta
      using VariableDelta = ::VariableDelta<Law,nDim,Basis>;

   // variables
      std::array<Types::Real,nVar<Law,nDim>> var{};

   // default, copy and move constructors
      VariableSet() = default;
      VariableSet( const VariableSet&  ) = default;
      VariableSet(       VariableSet&& ) = default;

   // only explicit construction from VariableDelta
      explicit VariableSet( const VariableDelta&  dq ) noexcept : var(dq.var) {}
      explicit VariableSet(       VariableDelta&& dq ) noexcept : var(std::move(dq.var)) {}

   // initialiser list constructor, must be a length nVar list of Types::Real
      template<typename... T>
         requires   Same<Types::Real,T...>
                 && is_pack_of_n<nVar<Law,nDim>,T...>::value
      VariableSet( T... r ) noexcept : var{r...} {}

   // copy/move assignment
      VariableSet& operator=( const VariableSet&  ) = default;
      VariableSet& operator=(       VariableSet&& ) = default;

   // accessors
            Types::Real& operator[]( const int i )       { return var[i]; }
      const Types::Real& operator[]( const int i ) const { return var[i]; }

   // in-place arithmetic
      VariableSet& operator+=( const VariableDelta& dq0 );
      VariableSet& operator-=( const VariableDelta& dq0 );
      VariableSet& operator =( const Types::Real a );
  };

/*
 * send each element of VariableSet to stream, seperated by a single space
 */
   template<LawType Law, int nDim, BasisType<Law> Basis>
   std::ostream& operator<<( std::ostream& os, const VariableSet<Law,nDim,Basis>& q );
                                              
/*
 * Vector of displacements for a hyperbolic conservation law (Law) in a particular basis for the phase space (Basis) with (nDim) spatial dimensions
 * VariableDelta is a displacement in an affine space, with VariableSet<Law,nDim,Basis> being points in this space
 */
   template<LawType Law, int nDim, BasisType<Law> Basis>
   struct VariableDelta
  {
   // equivalent VariableDelta
      using VariableSet = ::VariableSet<Law,nDim,Basis>;

   // variables
      std::array<Types::Real,nVar<Law,nDim>> var{};

   // default, copy and move constructors
      VariableDelta() = default;
      VariableDelta( const VariableDelta&  ) = default;
      VariableDelta(       VariableDelta&& ) = default;

   // only explicit conversion from Delta
      explicit VariableDelta( const VariableSet&  q ) noexcept : var(q.var) {};
      explicit VariableDelta(       VariableSet&& q ) noexcept : var(std::move(q.var)) {};

   // initialiser list constructor, must be a length nVar list of Types::Real
      template<typename... T>
         requires   Same<Types::Real,T...>
                 && is_pack_of_n<nVar<Law,nDim>,T...>::value
      VariableDelta( T... r ) noexcept : var{r...} {}

   // copy/move assignment
      VariableDelta& operator=( const VariableDelta&  ) = default;
      VariableDelta& operator=(       VariableDelta&& ) = default;

   // accessors
            Types::Real& operator[]( const int i )       { return var[i]; }
      const Types::Real& operator[]( const int i ) const { return var[i]; }

   // in-place arithmetic
      VariableDelta& operator+=( const VariableDelta& dq0 );
      VariableDelta& operator-=( const VariableDelta& dq0 );
      VariableDelta& operator*=( const Types::Real a );
      VariableDelta& operator/=( const Types::Real a );
      VariableDelta& operator =( const Types::Real a );
  };

/*
 * send each element of VariableDelta to stream, seperated by a single space
 */
   template<LawType Law, int nDim, BasisType<Law> Basis>
   std::ostream& operator<<( std::ostream& os, const VariableDelta<Law,nDim,Basis>& dq );
                                              
/*
 * A flux, and associated spectral radius in the phase space of conservation law Law in nDim spatial dimensions
 * flux is a VariableDelta in Conserved variables to ensure correct shock speeds
 * spectral radius is the largest eigenvalue of the flux jacobian, scaled by the area of the cell face the flux is over
 */
   template<LawType Law, int nDim>
   struct FluxResult
  {
      using FluxType = VariableDelta<Law,nDim,BasisType<Law>::Conserved>;

      FluxType      flux{};
      Types::Real lambda{};

   // default, copy and move constructors
      FluxResult() = default;
      FluxResult( const FluxResult&  ) = default;
      FluxResult(       FluxResult&& ) = default;

   // copy/move assignment
      FluxResult& operator=( const FluxResult&  ) = default;
      FluxResult& operator=(       FluxResult&& ) = default;

   // construct from flux and spectral radius
      FluxResult( const FluxType&  f, const Types::Real l ) noexcept : flux(f), lambda(l) {};
      FluxResult(       FluxType&& f, const Types::Real l ) noexcept : flux(std::move(f)), lambda(l) {};

   // in-place arithmetic
      // fluxes are added or subtracted
      // spectral radius is only added
      FluxResult& operator+=( const FluxResult& fr );
      FluxResult& operator-=( const FluxResult& fr );
  };


// ---------- exact physical flux ----------

/*
 * calculating exactFlux from a (minimally implemented) VariableSet just wraps calculating exactFlux from a State
 * this means that only calculation from a State needs implementing for each LawType
 */
   template<LawType Law, ImplementedVarSet VarT>
      requires   SameLaw<VarT,law_constant<Law>>
              && ImplementedLawType<Law>
   fluxresult_t<VarT> exactFlux( const Species<Law>& species, const Geometry::Direction<dim_of_v<VarT>>& normal, const VarT& q0 )
  {
      return exactFlux( species, normal, setToState( species, q0 ) );
  }

/*
 * Calculate the magnitude of velocity along given direction
 *    assumes State<Law,nDim> has a .velocity(int) member
 */
   template<LawType Law, int nDim>
   Types::Real projectedVelocity( const Geometry::Direction<nDim>& dir, const State<Law,nDim>& state )
  {
      Types::Real un=0;
      for( int i=0; i<nDim; i++ )
     {
         un+= dir[i]*state.velocity(i);
     }
      return un;
  }

/*
 * for each vector quantity in VarT, rotate the underlying spatial basis to the given metric
 *    uses assumption that vector quantities are always contiguous first elements of VariableSet
 */
   template<LawType Law, int nDim, ImplementedVarSet VarT>
      requires   SameLaw<VarT,law_constant<Law>>
              && SameDim<VarT,dim_constant<nDim>>
   VarT rotateToMetric( const Species<Law>& species, const Geometry::Metric<nDim>& metric, const VarT& q0 )
  {
      constexpr int nVec = nVectorQuantities<Law>;

      VarT q1(q0);
      for( int v=0; v<nVec; v++ )
     {
         const int off=v*nDim;
         for( int i=0; i<nDim; i++ )
        {
            q1[off+i]=0;
            for( int j=0; j<nDim; j++ )
           {
               q1[off+i] += metric[i][j]*q0[off+j];
           }
        }
     }
      return q1;
  }


// ---------- explicit transformation functions ----------

/*
 * transformations to the same type just return a copy
 */
   template<ImplementedVarSet VarT, LawType Law>
      requires SameLaw<VarT,law_constant<Law>>
   VarT set2Set( const Species<Law>& species, const VarT& q0 )
  {
      return VarT(q0);
  }

/*
 * transformations between two VariableSets of different bases (SrcT & DstT) defaults to transforming SrcT->State->DstT
 *    this means that only transformations to/from a State needs implementing for each (Law,Basis) pair
 */
   template<ImplementedVarSet DstT, ImplementedVarSet SrcT, LawType Law>
      requires   SameLaw<DstT,SrcT,law_constant<Law>>
              && SameDim<DstT,SrcT>

   DstT set2Set( const Species<Law>& species, const SrcT& q0 )
  {
      return state2Set<DstT>( species, set2State( species, q0 ) );
  }

/*
 * transformations between two VariableDeltas of different bases must be around a particular point in the phase space. If this point is defined by a VariableSet, it will be transformed to a State before transforming the VariableDeltas.
 *    this means that Delta transformations only need defining around a State for each Basis pair
 */
   template<ImplementedVarDelta DstT, ImplementedVarDelta SrcT, ImplementedVarSet VarT, LawType Law>
      requires   SameLaw<DstT,SrcT,VarT,law_constant<Law>>
              && SameDim<DstT,SrcT,VarT>
   DstT delta2Delta( const Species<Law>& species, const VarT& q0, const SrcT& dq0 )
  {
      return delta2Delta( species, set2State( species, q0 ), dq0 );
  }

// template<ImplementedVarDelta DstT, ImplementedVarDelta SrcT, LawType Law, int nDim>
//    requires   SameLaw<DstT,SrcT,law_constant<Law>>
//            && SameDim<DstT,SrcT,dim_constant<nDim>>
// DstT delta2Delta( const Species<Law>&, const State<Law,nDim>& state, SrcT& dq0 )
//{
//    return flux2Delta( species, state, delta2Flux( species, state, dq0 ) );
//}


// ---------- generic flux functions ----------

/*
 * Base class for CRTP to give the same interface for all flux functions.
 *    provides interface to call flux with any variable set basis. This automatically transforms to State
 *    this allows each flux function to be implemented only with State argument, but called with any VariableSet
 */
   template<typename Flux, LawType Law>
      requires   ImplementedLawType<Law>
 //           && FluxImplementation<Flux,Law>
   struct FluxInterface
  {
      template<ImplementedVarSet VarSetT, int nDim>
         requires   SameLaw<law_constant<Law>, VarSetT>
                 && SameDim<dim_constant<nDim>,VarSetT>
      FluxResult<Law,nDim> operator()( const Species<Law>&         species,
                                       const Geometry::Surface<nDim>& face,
                                       const VarSetT&                   ql,
                                       const VarSetT&                   qr ) const
     {
         return static_cast<const Flux*>(this)->flux( species, face,
                                                set2State( species, ql ),
                                                set2State( species, qr ) );
     }

      template<int nDim>
      FluxResult<Law,nDim> operator()( const Species<Law>&         species,
                                       const Geometry::Surface<nDim>& face,
                                       const State<Law,nDim>&           sl,
                                       const State<Law,nDim>&           sr ) const
     {
         return static_cast<const Flux*>(this)->flux( species, face, sl, sr );
     }
  };

/*
 * Interface flux using a central average of the left/right exact fluxes
 *    Uses CRTP to provide interface to call flux with any VariableSet basis
 */
   template<LawType Law>
      requires ImplementedLawType<Law>
   struct CentralFlux : FluxInterface<CentralFlux<Law>,Law>
  {
      template<int nDim>
      static FluxResult<Law,nDim> flux( const Species<Law>&         species,
                                        const Geometry::Surface<nDim>& face,
                                        const State<Law,nDim>&           sl,
                                        const State<Law,nDim>&           sr );
  };

/*
 * Interface flux using a central average of the exact fluxes with scalar dissipation proportional to the spectral radius of the central flux
 *    Uses CRTP to provide interface to call flux with any VariableSet basis
 */
   template<LawType Law>
      requires ImplementedLawType<Law>
   struct RusanovFlux : FluxInterface<RusanovFlux<Law>,Law>
  {
      template<int nDim>
      static FluxResult<Law,nDim> flux( const Species<Law>&         species,
                                        const Geometry::Surface<nDim>& face,
                                        const State<Law,nDim>&           sl,
                                        const State<Law,nDim>&           sr );
  };


// ---------- implementation files ----------

# include <conservationLaws/base/variables/variableSet.ipp>
# include <conservationLaws/base/variables/variableDelta.ipp>
# include <conservationLaws/base/variables/variableArithmetic.ipp>
# include <conservationLaws/base/variables/variableArrays.ipp>
# include <conservationLaws/base/fluxResult.ipp>

# include <conservationLaws/base/fluxes/centralFlux.ipp>
# include <conservationLaws/base/fluxes/rusanovFlux.ipp>


// ---------- example specialisations ----------

/*
// must define overload for each LawType
   template<int nDim>
   FluxResult<Euler,nDim> exactFlux( Species<Euler>, Geometry::Surface<nDim>, State<Euler,nDim> );

// must define pair of function overloads for transforming each set to/from a state
   template< EulerConservedVariables DstT, EulerState SrcT >
      requires SameDim<DstT,SrcT>
   DstT transform( Species<Euler>, SrcT );

   template< EulerState DstT, EulerConservedVariables SrcT >
      requires SameDim<DstT,SrcT>
   DstT transform( Species<Euler>, SrcT );

// must define pair of transformations between deltas at particular state
   template< EulerConservedDelta DstT, EulerPrimitiveDelta SrcT >
      requires SameDim<DstT,SrcT>
   DstT transform( Species<Euler>, State<Euler,SrcT::nDim>, SrcT, Geometry::Metric<SrcT::nDim> );

   template< EulerPrimitiveDelta DstT, EulerConservedDelta SrcT >
      requires SameDim<DstT,SrcT>
   DstT transform( Species<Euler>, State<Euler,SrcT::nDim>, SrcT, Geometry::Metric<SrcT::nDim> );

// can define function overloads for transforming directly between sets for improved performance
   template< EulerConservedVariables DstT, EulerPrimitiveVariables SrcT >
      requires SameDim<DstT,SrcT>
   DstT transform( Species<Euler>, SrcT );

   template< EulerConservedVariables DstT, EulerPrimitiveVariables SrcT >
      requires SameDim<DstT,SrcT>
   DstT transform( Species<Euler>, SrcT );
*/

