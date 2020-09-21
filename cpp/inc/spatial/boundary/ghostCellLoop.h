
# pragma once

# include <spatial/boundary/boundaryCondition.h>

# include <conservationLaws/base/base.h>
# include <solutionField/solutionField.h>

# include <geometry/geometry.h>
# include <mesh/mesh.h>

# include <parallalg/array.h>

# include <cassert>


/*
 * Accumulate cell residual contributions over ghost cell boundary faces for 1D domain
 */
   template<LawType                   Law,
            BoundaryType<Law>      BCType,
            ImplementedVarSet     SolVarT,
            ImplementedVarDelta   SolDelT,
            typename              FluxRes,
            typename           UpdateFunc,
            typename        HighOrderFlux,
            floating_point           Real>
      requires   ConsistentTypes<Law,
                                 1,
                                 Real,
                                 SolVarT,
                                 SolDelT>
              && std::is_same_v<FluxRes,
                                fluxresult_t<SolVarT>>
   void boundaryResidual( const BoundaryCondition<Law,BCType,UpdateFunc>&,
                          const size_t                     boundaryId,
                          const Mesh<1,Real>&                    mesh,
                          const HighOrderFlux&                 hoflux,
                          const Species<Law,Real>&            species,
                          const SolutionField<SolVarT,1>&           q,
                          const par::Array<std::array<SolDelT,1>,1>&          dq,
                                par::Array<   FluxRes,1>&         res )
  {
   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() == dq.shape() );
      assert( mesh.cells.shape() == q.interior.shape() );

   // boundary condition type selection was correct?
      assert( q.bcTypes[boundaryId] == BCType );

   // first order boundaries
      if( boundaryId==0 ) // left boundary
     {
         const par::Idx<1> ip{0};
         const par::Idx<1> ir{0};
         const par::Idx<1> ib{0};
         const FluxRes fr = hoflux( species,
                                    surface( mesh.nodes[ip] ),
                                    SolDelT{},         dq[ir][0],
//                                  SolDelT{},      SolDelT{},
                                    q.boundary[0][ib], q.interior[ir] );
         res[ir]+=fr;
     }
      else if( boundaryId==1 ) // right boundary
     {
         const size_t nc = q.interior.shape(0);

         const par::Idx<1> ip{nc};
         const par::Idx<1> il{nc-1};
         const par::Idx<1> ib{0};

         const FluxRes fr = hoflux( species,
                                    surface( mesh.nodes[ip] ),
                                    dq[il][0],      SolDelT{},
//                                  SolDelT{},      SolDelT{},
                                    q.interior[il], q.boundary[1][ib] );
         res[il]-=fr;
     }
      else{ assert( false && "invalid boundary id for 1D ghost cell flux, must be 0 or 1" ); }

      return;
  }

/*
 * Accumulate cell residual contributions over ghost cell boundary faces for 2D domain
 */
   template<LawType                   Law,
            BoundaryType<Law>      BCType,
            ImplementedVarSet     SolVarT,
            ImplementedVarDelta   SolDelT,
            typename              FluxRes,
            typename           UpdateFunc,
            typename        HighOrderFlux,
            floating_point           Real>
      requires   ConsistentTypes<Law,
                                 2,
                                 Real,
                                 SolVarT,
                                 SolDelT>
              && std::is_same_v<FluxRes,
                                fluxresult_t<SolVarT>>
   void boundaryResidual( const BoundaryCondition<Law,BCType,UpdateFunc>&,
                          const size_t                        boundaryId,
                          const Mesh<2,Real>&                       mesh,
                          const HighOrderFlux&                    hoflux,
                          const Species<Law,Real>&               species,
                          const SolutionField<SolVarT,2>&              q,
                          const par::Array<std::array<SolDelT,2>,2>&  dq,
                                par::Array<FluxRes,2>&               res )
  {
   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() == dq.shape() );
      assert( mesh.cells.shape() == q.interior.shape() );

   // boundary condition type selection was correct?
      assert( q.bcTypes[boundaryId] == BCType );

      const size_t ni = q.interior.shape(0);
      const size_t nj = q.interior.shape(1);
      const size_t bID= boundaryId;

      if( boundaryId==0 ) // left face
     {
         const size_t i=0;
         for( size_t j=0; j<nj; ++j )
        {
         // interior cell  index
            const par::Idx<2> ic{i,j};

         // boundary cell index
            const par::Idx<2> ib{j,0};
   
         // face node indices
            const par::Idx<2> ip0{i,j  };
            const par::Idx<2> ip1{i,j+1};
   
            const FluxRes fr = hoflux( species,
                                       surface( mesh.nodes[ip0],
                                                mesh.nodes[ip1] ),
                                       SolDelT{},           dq[ic][0],
                                       q.boundary[bID][ib], q.interior[ic] );
            res[ic]+=fr;
        }
     }
      else if( boundaryId==1 ) // right face
     {
         const size_t i=ni-1;
         for( size_t j=0; j<nj; ++j )
        {
         // interior cell  index
            const par::Idx<2> ic{i,j};

         // boundary cell index
            const par::Idx<2> ib{j,0};
   
         // face node indices
            const par::Idx<2> ip0{i+1,j  };
            const par::Idx<2> ip1{i+1,j+1};
   
            const FluxRes fr = hoflux( species,
                                       surface( mesh.nodes[ip0],
                                                mesh.nodes[ip1] ),
                                       dq[ic][0],      SolDelT{},       
                                       q.interior[ic], q.boundary[bID][ib] );
            res[ic]-=fr;
        }
     }
      else if ( boundaryId==2 ) // bottom face
     {
         const size_t j=0;
         for( size_t i=0; i<ni; ++i )
        {
         // interior cell  index
            const par::Idx<2> ic{i,j};

         // boundary cell index
            const par::Idx<2> ib{i,0};
   
         // face node indices
            const par::Idx<2> ip0{i+1,j};
            const par::Idx<2> ip1{i  ,j};
   
            const FluxRes fr = hoflux( species,
                                       surface( mesh.nodes[ip0],
                                                mesh.nodes[ip1] ),
                                       SolDelT{},           dq[ic][0],
                                       q.boundary[bID][ib], q.interior[ic] );
            res[ic]+=fr;
        }
     }
      else if ( boundaryId==3 ) // top face
     {
         const size_t j=nj-1;
         for( size_t i=0; i<ni; ++i )
        {
         // interior cell  index
            const par::Idx<2> ic{i,j};

         // boundary cell index
            const par::Idx<2> ib{i,0};
   
         // face node indices
            const par::Idx<2> ip0{i+1,j+1};
            const par::Idx<2> ip1{i  ,j+1};
   
            const FluxRes fr = hoflux( species,
                                       surface( mesh.nodes[ip0],
                                                mesh.nodes[ip1] ),
                                        dq[ic][0],      SolDelT{},        
                                        q.interior[ic], q.boundary[bID][ib] );
            res[ic]-=fr;
        }
     }
      else{ assert( false && "invalid boundary id for 2D ghost cell flux, must be 0,1,2 or 3" ); }

      return;
  }

