
# pragma once

# include <spatial/boundary/boundaryCondition.h>

# include <conservationLaws/base/base.h>
# include <solutionField/solutionField.h>

# include <geometry/geometry.h>
# include <mesh/mesh.h>

# include <parallalg/array.h>

# include <utils/concepts.h>

# include <cassert>


/*
 * Accumulate cell residual contributions over periodic boundary face for 1D domain
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
   void boundaryResidual( const BoundaryCondition<Law,BoundaryType<Law>::Periodic>&,
                          const size_t                           boundaryId,
                          const Mesh<1,Real>&                          mesh,
                          const HighOrderFlux&                       hoflux,
                          const Species<Law,Real>&                  species,
                          const SolutionField<SolVarT,1>&                 q,
                          const par::DualArray<std::array<SolDelT,1>,1>& dq,
                                par::DualArray<FluxRes,1>&              res )
  {
   // valid boundary?
      assert( (boundaryId==0 or boundaryId==1)
              && "invalid boundary id for 1D boundary condition, must be 0 or 1"  );

   // boundary condition type selection was correct?
      assert( q.bcTypes[boundaryId] == BoundaryType<Law>::Periodic );

   // flux contribution from periodic face should only be calculated from 'lower' face
      if( boundaryId==1 ){ return; }

   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() == dq.shape() );
      assert( mesh.cells.shape() == q.interior.shape() );

   // periodic boundary
      const size_t nc = q.interior.shape(0);

      using CellIdx = typename SolutionField<SolVarT,1>::VarField::IdxType;
      using NodeIdx = typename Mesh<1,Real>::NodeArray::IdxType;

      const NodeIdx ip{0};
      const CellIdx il{nc-1};
      const CellIdx ir{0};

      const FluxRes fr = hoflux( species,
                                 surface( mesh.nodes[ip] ),
                                 dq[il][0],      dq[ir][0],
                                 q.interior[il], q.interior[ir] );
      res[il]-=fr;
      res[ir]+=fr;

      return;
  }

/*
 * Accumulate cell residual contributions over periodic boundary faces for 2D domain
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
   void boundaryResidual( const BoundaryCondition<Law,BoundaryType<Law>::Periodic>&,
                          const size_t                            boundaryId,
                          const Mesh<2,Real>&                           mesh,
                          const HighOrderFlux&                        hoflux,
                          const Species<Law,Real>&                   species,
                          const SolutionField<SolVarT,2>&                  q,
                          const par::DualArray<std::array<SolDelT,2>,2>&  dq,
                                par::DualArray<FluxRes,2>&               res )
  {
//    return;
   // valid boundary?
      assert((   (boundaryId==0)
              or (boundaryId==1)
              or (boundaryId==2)
              or (boundaryId==3))
             && "invalid boundary id for 1D boundary condition, must be 0 or 1"  );

   // boundary condition type selection was correct?
      assert( q.bcTypes[boundaryId] == BoundaryType<Law>::Periodic );

   // flux contribution from periodic face should only be calculated from 'lower' face
      if( (boundaryId==1) or (boundaryId==3) ){ return; }

   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() == dq.shape() );
      assert( mesh.cells.shape() == q.interior.shape() );

      const size_t ni = q.interior.shape(0);
      const size_t nj = q.interior.shape(1);

      using CellIdx = typename SolutionField<SolVarT,2>::VarField::IdxType;
      using NodeIdx = typename Mesh<2,Real>::NodeArray::IdxType;

      if( boundaryId==0 ) // periodic boundary along +/- i direction
     {
         for( size_t j=0; j<nj; ++j )
        {
            const size_t il=ni-1;
            const size_t ir=0;
   
         // cell left/right indices
            const CellIdx icl{il,j};
            const CellIdx icr{ir,j};
   
         // face node indices
            const NodeIdx ip0{0,j  };
            const NodeIdx ip1{0,j+1};
   
            const FluxRes fr = hoflux( species,
                                       surface( mesh.nodes[ip0],
                                                mesh.nodes[ip1] ),
                                       dq[icl][0],      dq[icr][0],
                                       q.interior[icl], q.interior[icr] );
            res[icl]-=fr;
            res[icr]+=fr;
        }
     }
      else if ( boundaryId==2 ) // periodic boundary along +/- j direction
     { 
         for( size_t i=0; i<ni; ++i )
        {
            const size_t jl=nj-1;
            const size_t jr=0;
   
         // cell left/right indices
            const CellIdx icl{i,jl};
            const CellIdx icr{i,jr};
   
         // face node indices
            const NodeIdx ip0{i+1,0};
            const NodeIdx ip1{i  ,0};
   
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

