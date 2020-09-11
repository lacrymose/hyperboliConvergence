
# pragma once

# include <limiters/limiter.h>

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <mesh/mesh.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>
# include <vector>

# include <cassert>

/*
 * 2D Higher order flux evaluation with general form and ghost cells at boundaries
 */

/*
   template<LawType Law, ImplementedVarSet SolVarT, ImplementedVarDelta SolDelT, typename HighOrderFlux, floating_point Real>
      requires ConsistentTypes<Law,1,Real,SolVarT,SolDelT>
   void residualCalc( const Mesh<1,Real>&                        mesh,
                      const HighOrderFlux&                     hoflux,
                      const std::array<SolVarT,2>&                qbc,
                      const par::Array<SolVarT,1>&                  q,
                      const par::Array<SolDelT,1>&                 dq,
                            par::Array<FluxResult<Law,1,Real>,1>& res )
  {
      using FluxRes = FluxResult<Law,1,Real>;

   // check mesh sizes match

      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() ==  dq.shape() );
      assert( mesh.cells.shape() ==   q.shape() );
      assert( mesh.nodes.shape() == par::nodeDims_from_cellDims(mesh.cells.shape()) );

      const size_t nc = mesh.cells.shape(0);

      par::fill( res, FluxRes{} );

   // accumulate cell residual contributions from the flux across each face
      for( size_t i=0; i<nc-1; i++ )
     {
         const size_t ip = i;
         const size_t il = i;
         const size_t ir = i+1;

         const FluxRes fr = hoflux( surface( mesh.nodes[{ip}] ),
                                    dq[{il}],dq[{ir}],
                                     q[{il}], q[{ir}] );

         res[{il}]-=fr;
         res[{ir}]+=fr;
     }

   // first order boundaries
   // left boundary
     {
         const size_t ip=0;
         const size_t ir=0;
         const FluxRes fr = hoflux( surface( mesh.nodes[{ip}] ),
                                    SolDelT{},dq[{ir}],
                                    qbc[0],    q[{ir}] );
         res[{ir}]+=fr;
     }

   // right boundary
     {
         const size_t ip=nc-1;
         const size_t il=nc-1;
         const FluxRes fr = hoflux( surface( mesh.nodes[{ip}] ),
                                    dq[{il}],SolDelT{},
                                     q[{il}],   qbc[1] );
         res[{il}]-=fr;
     }

      return;
  }
*/


/*
 * 2D Higher order flux evaluation with general form and periodic boundaries
 */

   template<LawType Law, ImplementedVarSet SolVarT, ImplementedVarDelta SolDelT, typename HighOrderFlux, floating_point Real>
      requires ConsistentTypes<Law,2,Real,SolVarT,SolDelT>
   void residualCalc( const Mesh<2,Real>&                        mesh,
                      const HighOrderFlux&                     hoflux,
                      const par::Array<SolVarT,2>&                  q,
                      const par::Array<std::array<SolDelT,2>,2>&   dq,
                            par::Array<FluxResult<Law,2,Real>,2>& res )
  {
      using FluxRes = FluxResult<Law,2,Real>;

      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() ==  dq.shape() );
      assert( mesh.cells.shape() ==   q.shape() );
      assert( mesh.nodes.shape() == par::nodeDims_from_cellDims( mesh.cells.shape() ) );

      par::fill( res, FluxRes{} );

      const size_t ni = mesh.cells.shape(0);
      const size_t nj = mesh.cells.shape(1);

   // accumulate cell residual contributions from fluxes across i-normal faces
      for( size_t i=0; i<ni-1; i++ )
     {
         for( size_t j=0; j<nj; j++ )
        {
         // cell left/right indices
            const par::Idx<2> icl{i,  j};
            const par::Idx<2> icr{i+1,j};

         // face node indices
            const par::Idx<2> ip0{i+1,j  };
            const par::Idx<2> ip1{i+1,j+1};

            const FluxRes fr = hoflux( surface( mesh.nodes[ip0],
                                                mesh.nodes[ip1] ),
                                       dq[icl][0],dq[icr][0],
                                        q[icl],    q[icr] );
            res[icl]-=fr;
            res[icr]+=fr;
        }
     }

   // accumulate cell residual contributions from fluxes across j-normal faces
      for( size_t i=0; i<ni; i++ )
     {
         for( size_t j=0; j<nj-1; j++ )
        {
         // cell left/right indices
            const par::Idx<2> icl{i,j  };
            const par::Idx<2> icr{i,j+1};

         // face node indices
            const par::Idx<2> ip0{i+1,j+1};
            const par::Idx<2> ip1{i  ,j+1};

            const FluxRes fr = hoflux( surface( mesh.nodes[ip0],
                                                mesh.nodes[ip1] ),
                                       dq[icl][1],dq[icr][1],
                                        q[icl],    q[icr] );
            res[icl]-=fr;
            res[icr]+=fr;
        }
     }

   // periodic boundary along +/- i direction
      for( size_t j=0; j<nj; j++ )
     {
         const size_t il=ni-1;
         const size_t ir=0;

      // cell left/right indices
         const par::Idx<2> icl{il,j};
         const par::Idx<2> icr{ir,j};

      // face node indices
         const par::Idx<2> ip0{0,j  };
         const par::Idx<2> ip1{0,j+1};

         const FluxRes fr = hoflux( surface( mesh.nodes[ip0],
                                             mesh.nodes[ip1] ),
                                    dq[icl][0],dq[icr][0],
                                     q[icl],    q[icr] );
         res[icl]-=fr;
         res[icr]+=fr;
     }

   // periodic boundary along +/- j direction
      for( size_t i=0; i<ni; i++ )
     {
         const size_t jl=nj-1;
         const size_t jr=0;

      // cell left/right indices
         const par::Idx<2> icl{i,jl};
         const par::Idx<2> icr{i,jr};

      // face node indices
         const par::Idx<2> ip0{i+1,0};
         const par::Idx<2> ip1{i  ,0};

         const FluxRes fr = hoflux( surface( mesh.nodes[ip0],
                                             mesh.nodes[ip1] ),
                                    dq[icl][1],dq[icr][1],
                                     q[icl],    q[icr] );
         res[icl]-=fr;
         res[icr]+=fr;
     }

      return;
  }

/*
 * return value overload
 */
   template<ImplementedVarSet SolVarT, ImplementedVarDelta SolDelT, typename HighOrderFlux, floating_point Real>
      requires ConsistentTypes<law_of_v<SolVarT>,2,Real,SolVarT,SolDelT>
   par::Array<fluxresult_t<SolVarT>,2>
      residualCalc( const Mesh<2,Real>&                        mesh,
                    const HighOrderFlux&                     hoflux,
                    const par::Array<SolVarT,2>&                  q,
                    const par::Array<std::array<SolDelT,2>,2>&   dq )
  {
      assert( mesh.shape ==  q.shape() );
      assert( mesh.shape == dq.shape() );
      par::Array<fluxresult_t<SolVarT>,2> res(q.shape());

      residualCalc( mesh,hoflux,q,dq,res );
      return res;
  }
