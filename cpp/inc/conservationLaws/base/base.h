
# pragma once

# include <conservationLaws/base/concepts.h>
# include <conservationLaws/base/type_traits.h>
# include <conservationLaws/base/declarations.h>

# include <geometry/geometry.h>

# include <types.h>

# include <type_traits>
# include <ostream>

// ---------- integral values ----------

/*
 * unique identifier for each Conservation Law
*/
   enum LawType struct { NoLaw,
                         ArtificialCompressibility,
                         Burgers,
                         Euler,
                         IsothermalGas,
                         Maxwell,
                         ScalarAdvection,
                         ShallowWater,
                         TrafficFlow };

/*
 * total dimension of the phase space for each Conservation Law for a given number of spatial dimensions
 */
   template<LawType Law, int nDim>
   constexpr int nVar = nScalarQuantities<Law> + nDim*nVectorQuantities<Law>;


// ---------- types for variables in conservation law phase space ----------

/*
 * Vector of variables for a hyperbolic conservation law (Law) in a particular basis for the phase space (Basis) with (nDim) spatial dimensions
 * VariableSet is a point in an affine space, with VariableDelta<Law,nDim,Basis> being displacements in this space
 */
   template<LawType Law, int nDim, BasisType<Law> Basis>
   struct VariableSet
  {
   // equivalent VariableDelta
      using VariableDelta = VariableDelta<Law,nDim,Basis>;

   // variables
      std::array<Types::Real,nVar<Law,nDim>> var{0};

   // default, copy and move constructors
      VariableSet() = default;
      VariableSet( const VariableSet&  ) = default;
      VariableSet(       VariableSet&& ) = default;

   // only explicit conversion from VariableDelta
      explicit VariableSet( const VariableDelta&  dq ) noexcept;
      explicit VariableSet(       VariableDelta&& dq ) noexcept;

   // copy/move assignment
      VariableSet& operator=( const VariableSet&  ) = default;
      VariableSet& operator=(       VariableSet&& ) = default;

   // accessors
      inline       Types::Real& operator[]( const int i )       { return var[i]; }
      inline const Types::Real& operator[]( const int i ) const { return var[i]; }

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
      using VariableSet = VariableSet<Law,nDim,Basis>;

   // variables
      std::array<Types::Real,nVar<Law,nDim>> var{0};

   // default, copy and move constructors
      VariableDelta() = default;
      VariableDelta( const VariableDelta&  ) = default;
      VariableDelta(       VariableDelta&& ) = default;

   // only explicit conversion from Delta
      explicit VariableDelta( const VariableSet&  dq ) noexcept;
      explicit VariableDelta(       VariableSet&& dq ) noexcept;

   // copy/move assignment
      VariableDelta& operator=( const VariableDelta&  ) = default;
      VariableDelta& operator=(       VariableDelta&& ) = default;

   // accessors
      inline       Types::Real& operator[]( const int i )       { return var[i]; }
      inline const Types::Real& operator[]( const int i ) const { return var[i]; }

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

      FluxType      flux;
      Types::Real lambda;

   // default, copy and move constructors
      FluxResult() = default;
      FluxResult( const FluxResult&  ) = default
      FluxResult( const FluxResult&& ) = default

   // copy/move assignment
      FluxResult& operator=( const FluxResult&  ) = default;
      FluxResult& operator=( const FluxResult&& ) = default;

   // construct from flux and spectral radius
      FluxResult( FluxType f, Types::Real l ) : flux(f) lambda(l) {};
  };


// ---------- exact physical flux ----------

/*
 * calculating exactFlux from a (minimally implemented) VariableSet just wraps calculating exactFlux from a State
 * this means that only calculation from a State needs implementing for each LawType
 */
   template<LawType Law, ImplementedVarSet VarT>
      requires   SameLaw<VarT,Law>
              && ImplementedLawType<Law>
   fluxresult_t<VarT> exactFlux( Species<Law> species, Geometry::Surface<dim_of_v<VarT>> face, VarT q0 )
  {
      using StateT = state_t<SrcT>;
      return exactFlux( species, face, transform<StateT>( species, q0 ) );
  }

// ---------- transformation functions ----------

/*
 * transformations between two VariableSets of different bases (SrcT & DstT) defaults to transforming SrcT->State->DstT
 * this means that only transformations to/from a State needs implementing for each (Law,Basis) pair
 */
   template<ImplementedVarSet DstT, ImplementedVarSet SrcT, LawType Law>
      requires   SameLaw<DstT,SrcT,Law>
              && SameDim<DstT,SrcT>
              && ImplementedLawType<Law>
   DstT transform( Species<Law> species, SrcT q0 )
  {
      using StateT = state_t<SrcT>;
      return transform<DstT>( species, transform<StateT>( species, q0 ) );
  }

/*
 * transformations between two VariableDeltas of different bases must be around a particular point in the phase space. If this point is defined by a VariableSet, it will be transformed to a State before transforming the VariableDeltas.
 * this means that Delta transformations only need defining around a State for each Basis pair
 */
   template<ImplementedVarDelta DstT, ImplementedVarDelta SrcT, ImplementedVarSet VarT, LawType Law>
      requires   SameLaw<DstT,SrcT,VarT,Law>
              && SameDim<DstT,SrcT,VarT>
              && ImplementedLawType<Law>
   DstT transform( Species<Law> species, VarT q0, SrcT dq0, Geometry::Metric<dim_of<DstT> n )
  {
      using StateT = state_t<SrcT>;
      return transform<DstT>( species, transform<StateT>( species, q0 ), dq0, n );
  }


// ---------- generic flux functions ----------

/*
 * Base class for CRTP to give the same interface for all flux functions.
 *    provides interface to call flux with any variable set basis. This automatically transforms to State
 *    this allows each flux function to be implemented only with State argument, but called with any VariableSet
 */
   template<typename Flux, LawType Law>
      requires   ImplementedLawType<Law>
              && FluxImplementation<Flux,Law>
   struct FluxInterface
  {
      template<ImplementedVarSet VarSetT, int nDim>
         requires   SameLaw<Law,VarSetT>
                 && SameDim<nDim,VarSetT>
      FluxResult<Law,nDim> operator()( const Species<Law>          species,
                                       const Geometry::Surface<nDim>& face,
                                       const VarSetT                    ql,
                                       const VarSetT                    qr ) const
     {
         return static_cast<Flux*>(this)->flux( species, face,
                                                transform<State<Law,nDim>( species, ql ),
                                                transform<State<Law,nDim>( species, qr ) );
     };

      template<int nDim>
      FluxResult<Law,nDim> operator()( const Species<Law>          species,
                                       const Geometry::Surface<nDim>& face,
                                       const State<Law,nDim>            sl,
                                       const State<Law,nDim>            sr ) const
     {
         return static_cast<Flux*>(this)->flux( species, face, sl, sr );
     };
  };

/*
 * Interface flux using a central average of the left/right exact fluxes
 *    Uses CRTP to provide interface to call flux with any VariableSet basis
 */
   template<LawType Law>
      requires ImplementedLawType<Law>
   struct CentralFlux : FluxInterface<CentralFlux,Law>
  {
      template<int nDim>
      static FluxResult<Law,nDim> flux( const Species<Law>&         species,
                                        const Geometry::Surface<nDim>& face,
                                        const State<Law,nDim>&           sl,
                                        const State<Law,nDim>&           sr ) const;
  };

/*
 * Interface flux using a central average of the exact fluxes with scalar dissipation proportional to the spectral radius of the central flux
 *    Uses CRTP to provide interface to call flux with any VariableSet basis
 */
   template<LawType Law>
      requires ImplementedLawType<Law>
   struct RusanovFlux : FluxInterface<RusanovFlux,Law>
  {
   // types needed for upwind diffusion calculation
      using ConservedSet   = VariableSet<  Law,nDim,BasisType<Law>::Conserved>;
      using ConservedDelta = VariableDelta<Law,nDim,BasisType<Law>::Conserved>;

      template<int nDim>
      static FluxResult<Law,nDim> flux( const Species<Law>&         species,
                                        const Geometry::Surface<nDim>& face,
                                        const State<Law,nDim>&           sl,
                                        const State<Law,nDim>&           sr ) const;
  };


// ---------- implementation files ----------

# include <conservationLaws/base/variables/variableSet.ipp>
# include <conservationLaws/base/variables/variableDelta.ipp>
# include <conservationLaws/base/variables/variableArithmetic.ipp>
# include <conservationLaws/base/variables/variableArrays.ipp>

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

