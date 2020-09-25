
# pragma once

# include <spatial/boundary/periodicLoop.h>
# include <spatial/boundary/ghostCellLoop.h>
# include <spatial/boundary/boundaryFluxLoop.h>
# include <spatial/boundary/boundaryCondition.h>

# include <solutionField/solutionField.h>
# include <conservationLaws/base/base.h>

# include <mesh/mesh.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <utils/utils.h>
# include <utils/concepts.h>

# include <tuple>
# include <cassert>

/*
 * Accumulate cell residuals from fluxes over all cell faces
 */
   template<LawType                   Law,
            int                      nDim,
            size_t                      N,
            ImplementedVarSet     SolVarT,
            ImplementedVarDelta   SolDelT,
            typename              FluxRes,
            typename        HighOrderFlux,
            typename...     BoundaryConds,
            floating_point           Real>
      requires   ConsistentTypes<Law,
                                 nDim,
                                 Real,
                                 SolVarT,
                                 SolDelT>
              && std::is_same_v<FluxRes,
                                fluxresult_t<SolVarT>>
              && N==nDim
   void residualCalc( const Mesh<nDim,Real>&                      mesh,
                      const HighOrderFlux&                      hoflux,
                      const std::tuple<BoundaryConds...>           bcs,
                      const Species<Law,Real>&                 species,
                      const SolutionField<SolVarT,nDim>&             q,
                      const par::Array<std::array<SolDelT,N>,nDim>& dq,
                            par::Array<   FluxRes,nDim>&           res )

  {
   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() ==  dq.shape() );
      assert( mesh.cells.shape() ==  q.interior.shape() );

      par::fill( res, FluxRes{} );

   // accumulate cell residual contributions from interior faces
      interiorResidual( mesh, hoflux, species, q,dq, res );

   // accumulate cell residual contributions from boundary faces
      boundaryResidual( mesh, hoflux, bcs, species, q,dq, res );

      return;
  }

/*
 * Accumulate cell residual contributions from interior faces in 1D domain
 */
   template<LawType                   Law,
            ImplementedVarSet     SolVarT,
            ImplementedVarDelta   SolDelT,
            typename              FluxRes,
            typename        HighOrderFlux,
            floating_point           Real>
      requires   ConsistentTypes<Law,
                                 1,
                                 Real,
                                 SolVarT,
                                 SolDelT>
              && std::is_same_v<FluxRes,
                                fluxresult_t<SolVarT>>
   void interiorResidual( const Mesh<1,Real>&               mesh,
                          const HighOrderFlux&            hoflux,
                          const Species<Law,Real>&       species,
                          const SolutionField<SolVarT,1>&      q,
                          const par::Array<std::array<SolDelT,1>,1>&     dq,
                                par::Array<FluxRes,1>&    res )
  {
   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() == dq.shape() );
      assert( mesh.cells.shape() == q.interior.shape() );

      const size_t nc = mesh.cells.shape(0);

   // accumulate cell residual contributions from the flux across each face
      for( size_t i=0; i<nc-1; ++i )
     {
         const par::Idx<1> ip{i};
         const par::Idx<1> il{i};
         const par::Idx<1> ir{i+1};

         const FluxRes fr = hoflux( species,
                                    surface( mesh.nodes[ip] ),
                                    dq[il][0],      dq[ir][0],
                                    q.interior[il], q.interior[ir] );
         res[il]-=fr;
         res[ir]+=fr;
     }

      return;
  }

/*
 * Accumulate cell residual contributions from interior faces in 2D domain
 */
   template<LawType                   Law,
            ImplementedVarSet     SolVarT,
            ImplementedVarDelta   SolDelT,
            typename              FluxRes,
            typename        HighOrderFlux,
            floating_point           Real>
      requires   ConsistentTypes<Law,
                                 2,
                                 Real,
                                 SolVarT,
                                 SolDelT>
              && std::is_same_v<FluxRes,
                                fluxresult_t<SolVarT>>
   void interiorResidual( const Mesh<2,Real>&                       mesh,
                          const HighOrderFlux&                    hoflux,
                          const Species<Law,Real>&               species,
                          const SolutionField<SolVarT,2>&              q,
                          const par::Array<std::array<SolDelT,2>,2>&  dq,
                                par::Array<FluxRes,2>&               res )
  {
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() == dq.shape() );
      assert( mesh.cells.shape() == q.interior.shape() );

      const size_t ni = mesh.cells.shape(0);
      const size_t nj = mesh.cells.shape(1);

   // accumulate cell residual contributions from fluxes across i-normal faces
# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t j=0; j<nj; ++j )
     {
         for( size_t i=0; i<ni-1; ++i )
        {
         // cell left/right indices
            const par::Idx<2> icl{i,  j};
            const par::Idx<2> icr{i+1,j};

         // face node indices
            const par::Idx<2> ip0{i+1,j  };
            const par::Idx<2> ip1{i+1,j+1};

            const FluxRes fr = hoflux( species,
                                       surface( mesh.nodes[ip0],
                                                mesh.nodes[ip1] ),
                                       dq[icl][0],      dq[icr][0],
                                       q.interior[icl], q.interior[icr] );
            res[icl]-=fr;
            res[icr]+=fr;
        }
     }

   // accumulate cell residual contributions from fluxes across j-normal faces
# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj-1; ++j )
        {
         // cell left/right indices
            const par::Idx<2> icl{i,j  };
            const par::Idx<2> icr{i,j+1};

         // face node indices
            const par::Idx<2> ip0{i+1,j+1};
            const par::Idx<2> ip1{i  ,j+1};

            const FluxRes fr = hoflux( species,
                                       surface( mesh.nodes[ip0],
                                                mesh.nodes[ip1] ),
                                       dq[icl][1],      dq[icr][1],
                                       q.interior[icl], q.interior[icr] );
            res[icl]-=fr;
            res[icr]+=fr;
        }
     }

      return;
  }

/*
 * Accumulate cell residuals from fluxes over boundary cell faces
 */
   template<LawType                   Law,
            int                      nDim,
            size_t                      N,
            ImplementedVarSet     SolVarT,
            ImplementedVarDelta   SolDelT,
            typename              FluxRes,
            typename        HighOrderFlux,
            typename...     BoundaryConds,
            floating_point           Real>
      requires   ConsistentTypes<Law,
                                 nDim,
                                 Real,
                                 SolVarT,
                                 SolDelT>
              && std::is_same_v<FluxRes,
                                fluxresult_t<SolVarT>>
              && N==nDim
   void boundaryResidual( const Mesh<nDim,Real>&                      mesh,
                          const HighOrderFlux&                      hoflux,
                          const std::tuple<BoundaryConds...>           bcs,
                          const Species<Law,Real>&                 species,
                          const SolutionField<SolVarT,nDim>&             q,
                          const par::Array<std::array<SolDelT,N>,nDim>& dq,
                                par::Array<   FluxRes,nDim>&           res )
  {
   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() == dq.shape() );
      assert( mesh.cells.shape() == q.interior.shape() );

   // if boundary condition type matches type of bc in tuple, calculate boundary residual
      const auto call_bc_resid = []( auto&&... args )
     {
         boundaryResidual( std::forward<decltype(args)>(args)... );
     };

   // for each boundary, calculate the residual from the matching boundary condition in the bc tuple
      for( unsigned int boundaryId=0; boundaryId<q.nBoundaries; ++boundaryId )
     {
         selectBoundaryCondition( q.bcTypes[boundaryId], bcs, call_bc_resid,
                                  boundaryId, mesh, hoflux, species, q, dq, res );
     }
      return;
  }

