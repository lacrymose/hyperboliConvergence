
# pragma once

# include <limiters/limiter.h>

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <vector>

# include <cassert>

/*
 * Higher order flux evaluation with general form and ghost cells at boundaries
 */

   template<LawType Law, ImplementedVarSet SolVarT, ImplementedVarDelta SolDelT, typename HighOrderFlux, floating_point Real>
      requires ConsistentTypes<Law,1,Real,SolVarT,SolDelT>
   void residualCalc( const geom::Mesh<1,Real>&                  mesh,
                      const HighOrderFlux                      hoflux,
                      const std::array<SolVarT,2>                 qbc,
                      const par::Array<SolVarT,1>&                  q,
                      const par::Array<SolDelT,1>&                 dq,
                            par::Array<FluxResult<Law,1,Real>,1>& res )
  {
      using FluxRes = FluxResult<Law,1,Real>;

   // check mesh sizes match

      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() ==  dq.shape() );
      assert( mesh.cells.shape() ==   q.shape() );
      assert( mesh.nodes.shape() == par::Shape<1>{mesh.cells.shape(0)+1} );

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


/*
 * Higher order flux evaluation with general form and periodic boundaries
 */

   template<LawType Law, ImplementedVarSet SolVarT, ImplementedVarDelta SolDelT, typename HighOrderFlux, floating_point Real>
      requires ConsistentTypes<Law,1,Real,SolVarT,SolDelT>
   void residualCalc( const geom::Mesh<1,Real>&                  mesh,
                      const HighOrderFlux                      hoflux,
                      const par::Array<SolVarT,1>&                  q,
                      const par::Array<SolDelT,1>&                 dq,
                            par::Array<FluxResult<Law,1,Real>,1>& res )
  {
      using FluxRes = FluxResult<Law,1,Real>;

   // check mesh sizes match

      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() ==  dq.shape() );
      assert( mesh.cells.shape() ==   q.shape() );
      assert( mesh.nodes.shape() == par::Shape<1>{mesh.cells.shape(0)+1} );

      const size_t nc = mesh.cells.shape(0);

      par::fill( res, FluxRes{} );

   // accumulate cell residual contributions from the flux across each face
      for( size_t i=0; i<nc-1; i++ )
     {
         const size_t ip=i;
         const size_t il=i;
         const size_t ir=i+1;

         const FluxRes fr = hoflux( surface( mesh.nodes[{ip}] ),
                                    dq[{il}],dq[{ir}],
                                     q[{il}], q[{ir}] );

         res[{il}]-=fr;
         res[{ir}]+=fr;
     }

   // periodic boundary
     {
         const size_t ip=0;
         const size_t il=nc-1;
         const size_t ir=0;

         const FluxRes fr = hoflux( surface( mesh.nodes[{ip}] ),
                                    dq[{il}],dq[{ir}],
                                     q[{il}], q[{ir}] );
         res[{il}]-=fr;
         res[{ir}]+=fr;
     }
      return;
  }

