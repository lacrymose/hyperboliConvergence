
# pragma once

# include <conservationLaws/base/concepts.h>
# include <conservationLaws/base/type-traits.h>
# include <conservationLaws/base/declarations.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>

# include <type_traits>
# include <ostream>


// ---------- types for variables in conservation law phase space ----------

/*
 * Vector of variables for a hyperbolic conservation law (Law) in a particular basis for the phase space (Basis) with (nDim) spatial dimensions
 * VariableSet is a point in an affine space, with VariableDelta<Law,nDim,Basis> being displacements in this space
 */
   template<LawType Law, int nDim, BasisType<Law> Basis, floating_point Real>
   struct VariableSet : AffinePointBase<nVar<Law,nDim>,
                                        VariableSet<  Law,nDim,Basis,Real>,
                                        VariableDelta<Law,nDim,Basis,Real>,
                                        Real>
  {
      using AffinePointBase<nVar<Law,nDim>,
                            VariableSet<  Law,nDim,Basis,Real>,
                            VariableDelta<Law,nDim,Basis,Real>,
                            Real>::AffinePointBase;
  };


/*
 * Vector of displacements for a hyperbolic conservation law (Law) in a particular basis for the phase space (Basis) with (nDim) spatial dimensions
 * VariableDelta is a displacement in an affine space, with VariableSet<Law,nDim,Basis> being points in this space
 */
   template<LawType Law, int nDim, BasisType<Law> Basis, floating_point Real>
   struct VariableDelta : AffineDeltaBase<nVar<Law,nDim>,
                                          VariableSet<  Law,nDim,Basis,Real>,
                                          VariableDelta<Law,nDim,Basis,Real>,
                                          Real>
  {
      using AffineDeltaBase<nVar<Law,nDim>,
                            VariableSet<  Law,nDim,Basis,Real>,
                            VariableDelta<Law,nDim,Basis,Real>,
                            Real>::AffineDeltaBase;
  };


/*
 * A flux, and associated spectral radius in the phase space of conservation law Law in nDim spatial dimensions
 * flux is a VariableDelta in Conserved variables to ensure correct shock speeds
 * spectral radius is the largest eigenvalue of the flux jacobian, scaled by the area of the cell face the flux is over
 */
   template<LawType Law, int nDim, floating_point Real>
   struct FluxResult
  {
      using FluxType = VariableDelta<Law,nDim,BasisType<Law>::Conserved,Real>;

      FluxType   flux;
      Real     lambda;

   // default, copy and move constructors
      FluxResult() = default;
      FluxResult( const FluxResult&  ) = default;
      FluxResult(       FluxResult&& ) = default;

   // copy/move assignment
      FluxResult& operator=( const FluxResult&  ) = default;
      FluxResult& operator=(       FluxResult&& ) = default;

   // construct from flux and spectral radius
      FluxResult( const FluxType&  f, const Real l ) noexcept : flux(f), lambda(l) {};
      FluxResult(       FluxType&& f, const Real l ) noexcept : flux(std::move(f)), lambda(l) {};

   // in-place arithmetic
      // fluxes are added or subtracted
      // spectral radius is only added
      FluxResult& operator+=( const FluxResult& fr );
      FluxResult& operator-=( const FluxResult& fr );
  };

/*
 * The wavespeeds (eigenvalues) associated with each characteristic field (eigenvectors)
 *    used for flux difference splitting (eg roe) schemes and preconditioning
 *    can be multiplied by a characteristic VariableDelta to convect
 */
   template<LawType Law, int nDim, floating_point Real>
   struct WaveSpeeds
  {
      std::array<Real,nVar<Law,nDim>> speeds;

   // accessors
            Real& operator[]( const int i )       { return speeds[i]; }
      const Real& operator[]( const int i ) const { return speeds[i]; }
  };

/*
 * Calculate wavespeeds for 1D split flux-jacobian aligned with face
 */
   template<LawType Law, int nDim, floating_point Real>
   WaveSpeeds<Law,nDim,Real> wavespeeds( const Species<Law,Real>&     species,
                                         const geom::Surface<nDim,Real>& face,
                                         const State<Law,nDim,Real>&    state )
  {
      return wavespeeds( species, face.metric[0], state );
  }


/*
 * multiply the jump in each characteristic field by the corresponding wavespeed
 */
   template<LawType Law, int nDim, floating_point Real>
   VariableDelta<Law,nDim,BasisType<Law>::Characteristic,Real>
      operator*( const WaveSpeeds<   Law,nDim,Real>&                               speeds,
                 const VariableDelta<Law,nDim,BasisType<Law>::Characteristic,Real>& waves )
  {
      constexpr int nv = nVar<Law,nDim>;
      VariableDelta<Law,nDim,BasisType<Law>::Characteristic,Real> result;
      for( int i=0; i<nv; ++i ){ result[i]=speeds[i]*waves[i]; }
      return result;
  }


// ---------- exact physical flux ----------

/*
 * calculating exactFlux from a (minimally implemented) VariableSet just wraps calculating exactFlux from a State
 * this means that only calculation from a State needs implementing for each LawType
 */
   template<LawType Law, ImplementedVarSet VarT, floating_point Real>
      requires   SameLaw<VarT,law_constant<Law>>
              && SameFPType<VarT,Real>
              && ImplementedLawType<Law>
   fluxresult_t<VarT> exactFlux( const Species<Law,Real>&                   species,
                                 const geom::Direction<dim_of_v<VarT>,Real>& normal,
                                 const VarT&                                     q0 )
  {
      return exactFlux( species, normal, setToState( species, q0 ) );
  }

/*
 * Calculate the magnitude of velocity along given direction
 *    assumes State<Law,nDim,Real> has a .velocity(int) member
 */
   template<LawType Law, int nDim, floating_point Real>
   Real projectedVelocity( const geom::Direction<nDim,Real>& dir,
                           const State<Law,nDim,Real>&     state )
  {
      Real un=0;
      for( int i=0; i<nDim; i++ )
     {
         un+= dir[i]*state.velocity(i);
     }
      return un;
  }

   template<LawType Law, int nDim, floating_point Real>
   Real projectedVelocity( const geom::Surface<nDim,Real>& face,
                           const State<Law,nDim,Real>&    state )
  {
      return projectedVelocity( face.metric[0], state );
  }

/*
 * for each vector quantity in VarT, rotate the underlying spatial basis from the background (identity) metric to the given metric
 *    assumes that vector quantities are always contiguous first elements of VariableSet
 */
   template<int nDim, typename VarT, floating_point Real>
      requires   (   ImplementedVarSet<  VarT>
                  || ImplementedVarDelta<VarT> )
              && ConsistentTypes<law_of_v<VarT>,nDim,Real,VarT>
   VarT rotateToMetric( const geom::Metric<nDim,Real>& metric,
                        const VarT&                        q0 )
  {
      constexpr LawType Law = law_of_v<VarT>;
      constexpr int    nVec = nVectorQuantities<Law>;

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

/*
 * for each vector quantity in VarT, rotate the underlying spatial basis from the given metric to the background (identity) metric
 *    assumes that vector quantities are always contiguous first elements of VariableSet
 */
   template<int nDim, typename VarT, floating_point Real>
      requires   (   ImplementedVarSet<  VarT>
                  || ImplementedVarDelta<VarT> )
              && ConsistentTypes<law_of_v<VarT>,nDim,Real,VarT>
   VarT rotateFromMetric( const geom::Metric<nDim,Real>& metric,
                          const VarT&                        q0 )
  {
      return rotateToMetric( geom::transpose( metric ), q0 );
  }

/*
 * overloads for rotating to a surface aligned metric
 */
   template<int nDim, typename VarT, floating_point Real>
      requires   (   ImplementedVarSet<  VarT>
                  || ImplementedVarDelta<VarT> )
              && ConsistentTypes<law_of_v<VarT>,nDim,Real,VarT>
   VarT rotateToFace( const geom::Surface<nDim,Real>& face,
                      const VarT&                       q0 )
  {
      return rotateToMetric( face.metric, q0 );
  }

   template<int nDim, typename VarT, floating_point Real>
      requires   (   ImplementedVarSet<  VarT>
                  || ImplementedVarDelta<VarT> )
              && ConsistentTypes<law_of_v<VarT>,nDim,Real,VarT>
   VarT rotateFromFace( const geom::Surface<nDim,Real>& face,
                        const VarT&                       q0 )
  {
      return rotateFromMetric( face.metric, q0 );
  }

/*
 * wrapper for spectral radius over a surface
 */
   template<LawType Law, int nDim, floating_point Real>
   Real spectralRadius( const geom::Surface<nDim,Real>& face,
                        const State<Law,nDim,Real>&     state )
  {
     return spectralRadius( face.metric[0], state );
  }

// ---------- explicit transformation functions ----------

/*
 * transformations between two VariableSets of different bases (SrcT & DstT) defaults to transforming SrcT->State->DstT
 *    this means that only transformations to/from a State needs implementing for each (Law,Basis) pair
 */
   template<ImplementedVarSet DstT,
            ImplementedVarSet SrcT,
            LawType            Law,
            floating_point    Real>
      requires   ConsistentTypes<Law,dim_of_v<SrcT>,Real,
                                 DstT,SrcT>
              && !(std::is_same_v<DstT,SrcT>)
   DstT set2Set( const Species<Law,Real>& species,
                 const SrcT&                   q0 )
  {
      return state2Set<DstT>( species, set2State( species, q0 ) );
  }

/*
 * transformations to the same type just return a copy
 */
   template<ImplementedVarSet VarT, LawType Law, floating_point Real>
      requires   SameLaw<VarT,law_constant<Law>>
              && SameFPType<VarT,Real>
   VarT set2Set( const Species<Law,Real>& species,
                 const VarT&                   q0 )
  {
      return VarT(q0);
  }


/*
 * transformations between two VariableDeltas of different bases (SrcT & DstT) defaults to transforming SrcT->ConservedDelta->DstT
 * transformations between two VariableDeltas at a particular point in the phase space, defined by a state
 *    this means that Delta transformations only need defining around a State for each Basis pair
 */
   template<ImplementedVarDelta DstT,
            ImplementedVarDelta SrcT,
            typename          StateT,
            LawType              Law,
            floating_point      Real>
      requires   is_State_v<StateT>
              && ConsistentTypes<Law,dim_of_v<SrcT>,Real,
                                 DstT,SrcT,StateT>
              && !(std::is_same_v<DstT,SrcT>)
   DstT delta2Delta( const Species<Law,Real>&      species,
                     const StateT&                   state,
                     const SrcT&                       dq0 )
  {
      constexpr int nDim = dim_of_v<SrcT>;
      using ConsDelT = VariableDelta<Law,nDim,BasisType<Law>::Conserved,Real>;

      return delta2Delta<DstT>( species,
                                state,
                                delta2Delta<ConsDelT>( species,
                                                       state,
                                                       dq0 ) );
  }

/*
 * If point in phase space is defined by a VariableSet, it will be transformed to a State before transforming the VariableDeltas.
*/
   template<ImplementedVarDelta DstT,
            ImplementedVarDelta SrcT,
            ImplementedVarSet   VarT,
            LawType              Law,
            floating_point      Real>
      requires   ConsistentTypes<Law,dim_of_v<SrcT>,Real,
                                 DstT,SrcT,VarT>
              && !(std::is_same_v<DstT,SrcT>)
   DstT delta2Delta( const Species<Law,Real>&      species,
                     const VarT&                        q0,
                     const SrcT&                       dq0 )
  {
      return delta2Delta<DstT>( species,
                                set2State( species, q0 ),
                                dq0 );
  }

/*
 * transformations to the same type just return a copy
 */
   template<ImplementedVarDelta DelT,
            typename            VarT,
            LawType              Law,
            floating_point      Real>
      requires   (   ImplementedVarSet<VarT>
                  || is_State_v<VarT> )
              && ConsistentTypes<Law,
                                 dim_of_v<DelT>,
                                 Real,
                                 DelT,
                                 VarT>
   DelT delta2Delta( const Species<Law,Real>&      species,
                     const VarT&                        q0,
                     const DelT&                       dq0 )
  {
      return DelT(dq0);
  }


   template<LawType            Law,
            ImplementedVarSet VarT,
            floating_point    Real>
      requires ConsistentTypes<Law,
                               dim_of_v<VarT>,
                               Real,
                               VarT>
   state_t<VarT> roeAverage( const Species<Law,Real>& species,
                             const VarT&                   ql,
                             const VarT&                   qr )
  {
      return roeAverage( species,
                         set2State( species, ql ),
                         set2State( species, qr ) );
  }


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
      template<ImplementedVarSet VarSetT, int nDim, floating_point Real>
         requires   SameLaw<law_constant<Law>, VarSetT>
                 && SameDim<dim_constant<nDim>,VarSetT>
              && SameFPType<VarSetT,Real>
      FluxResult<Law,nDim,Real> operator()( const Species<Law,Real>&     species,
                                            const geom::Surface<nDim,Real>& face,
                                            const VarSetT&                    ql,
                                            const VarSetT&                    qr ) const
     {
         return static_cast<const Flux*>(this)->flux( species, face,
                                                      set2State( species, ql ),
                                                      set2State( species, qr ) );
     }

      template<int nDim, floating_point Real>
      FluxResult<Law,nDim,Real> operator()( const Species<Law,Real>&     species,
                                            const geom::Surface<nDim,Real>& face,
                                            const State<Law,nDim,Real>&       sl,
                                            const State<Law,nDim,Real>&       sr ) const
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
      template<int nDim, floating_point Real>
      static FluxResult<Law,nDim,Real> flux( const Species<Law,Real>&     species,
                                             const geom::Surface<nDim,Real>& face,
                                             const State<Law,nDim,Real>&       sl,
                                             const State<Law,nDim,Real>&       sr );
  };

/*
 * Interface flux using a central average of the exact fluxes with scalar upwind dissipation proportional to the spectral radius of the central flux
 *    Uses CRTP to provide interface to call flux with any VariableSet basis
 */
   template<LawType Law>
      requires ImplementedLawType<Law>
   struct RusanovFlux : FluxInterface<RusanovFlux<Law>,Law>
  {
      template<int nDim, floating_point Real>
      static FluxResult<Law,nDim,Real> flux( const Species<Law,Real>&     species,
                                             const geom::Surface<nDim,Real>& face,
                                             const State<Law,nDim,Real>&       sl,
                                             const State<Law,nDim,Real>&       sr );
  };

/*
 * Interface dissipative flux of scalar dissipation proportional to the spectral radius of the central flux
 *    Uses CRTP to provide interface to call flux with any VariableSet basis
 */
   template<LawType Law>
      requires ImplementedLawType<Law>
   struct RusanovDissipation : FluxInterface<RusanovDissipation<Law>,Law>
  {
      template<int nDim, floating_point Real>
      static FluxResult<Law,nDim,Real> flux( const Species<Law,Real>&     species,
                                             const geom::Surface<nDim,Real>& face,
                                             const State<Law,nDim,Real>&       sl,
                                             const State<Law,nDim,Real>&       sr );
  };

/*
 * Interface flux using a central average of the exact fluxes with scalar upwind dissipation proportional to the spectral radius of the central flux
 *    Uses CRTP to provide interface to call flux with any VariableSet basis
 */
   template<LawType Law>
      requires ImplementedLawType<Law>
   struct RoeFlux : FluxInterface<RoeFlux<Law>,Law>
  {
      template<int nDim, floating_point Real>
      static FluxResult<Law,nDim,Real> flux( const Species<Law,Real>&     species,
                                             const geom::Surface<nDim,Real>& face,
                                             const State<Law,nDim,Real>&       sl,
                                             const State<Law,nDim,Real>&       sr );
  };



// ---------- implementation files ----------

# include <conservationLaws/base/fluxResult.ipp>

# include <conservationLaws/base/fluxes/centralFlux.ipp>
# include <conservationLaws/base/fluxes/rusanovFlux.ipp>
# include <conservationLaws/base/fluxes/roeFlux.ipp>

