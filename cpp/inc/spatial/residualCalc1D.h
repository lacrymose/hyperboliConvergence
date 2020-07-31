
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
 * Accumulate cell residual contributions from interior faces
 *    avoid interior face of cell on characteristic boundary (characteristic bc contributes all df/dn)
 */
   template<LawType                   Law,
            ImplementedVarSet     SolVarT,
            ImplementedVarDelta   SolDelT,
            typename              FluxRes,
            typename        HighOrderFlux,
            floating_point           Real>
      requires   ConsistentTypes<Law,1,Real,SolVarT,SolDelT>
              && std::is_same_v<FluxRes,fluxresult_t<SolVarT>>
   void interiorResidualCalc( const Mesh<1,Real>&                        mesh,
                              const HighOrderFlux&                     hoflux,
                              const std::array<BoundaryType<Law>,2>   bctypes,  // which bc type is each boundary?
                              const par::Array<SolVarT,1>&                  q,  // current solution
                              const par::Array<SolDelT,1>&                 dq,  // gradients
                                    par::Array<FluxRes,1>&                res )
  {
   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() ==  dq.shape() );
      assert( mesh.cells.shape() ==   q.shape() );

   // if characteristic boundary, avoid face first interior face
      const iLo = [&]()
     {
         if( bctypes[0]==BoundaryType<Law>::Characteristic ) { return 1; }
         else { return 0; }
     }();

      const iHi = [&]()
     {
         if( bctypes[0]==BoundaryType<Law>::Characteristic ) { return q.shape(0)-2; }
         else { return q.shape(0)-1; }
     }();

   // accumulate cell residual contributions from the flux across each face
      for( size_t i=iLo; i<iHi; ++i )
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
      return;
  }

/*
 * Accumulate cell residual contributions over periodic boundary faces
 */
   template<LawType                   Law,
            ImplementedVarSet     SolVarT,
            ImplementedVarDelta   SolDelT,
            typename              FluxRes,
            typename        HighOrderFlux,
            floating_point           Real>
      requires   ConsistentTypes<Law,1,Real,SolVarT,SolDelT>
              && std::is_same_v<FluxRes,fluxresult_t<SolVarT>>
   void periodicResidualCalc( const Mesh<1,Real>&                    mesh,
                              const HighOrderFlux&                 hoflux,
                              const size_t                     boundaryId,  // which bc type is each boundary?
                              const par::Array<SolVarT,1>&              q,  // current solution
                              const par::Array<SolDelT,1>&             dq,  // gradients
                              const par::Array<SolVarT,1>&             qb,
                                    par::Array<FluxRes,1>&            res )
  {
   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() ==  dq.shape() );
      assert( mesh.cells.shape() ==   q.shape() );

   // flux contribution from periodic face should only be calculated from 'lower' face
      assert( boundaryId==0 );

   // periodic boundary
      const size_t ip=0;
      const size_t il=nc-1;
      const size_t ir=0;

      const FluxRes fr = hoflux( surface( mesh.nodes[{ip}] ),
                                 dq[{il}],dq[{ir}],
                                  q[{il}], q[{ir}] );
      res[{il}]-=fr;
      res[{ir}]+=fr;

      return;
  }

/*
 * Accumulate cell residual contributions over faces with characteristic boundary condition
 */
   template<LawType                   Law,
            ImplementedVarSet     SolVarT,
            ImplementedVarDelta   SolDelT,
            typename              FluxRes,
            floating_point           Real>
      requires   ConsistentTypes<Law,1,Real,SolVarT,SolDelT>
              && std::is_same_v<FluxRes,fluxresult_t<SolVarT>>
   void characteristicResidualCalc( const Mesh<1,Real>&                    mesh,
                                    const size_t                     boundaryId,  // which bc type is each boundary?
                                    const par::Array<SolVarT,1>&              q,  // current solution
                                    const par::Array<SolDelT,1>&             dq,  // gradients
                                    const par::Array<SolVarT,1>&             qb,
                                          par::Array<FluxRes,1>&            res )
  {
   // check mesh sizes match
      assert( mesh.cells.shape() == res.shape() );
      assert( mesh.cells.shape() ==  dq.shape() );
      assert( mesh.cells.shape() ==   q.shape() );

   // characteristic boundary
      if( boundaryId==0 )
     {
         const par::Idx<1> ip{0};
         const par::Idx<1> ic{0};
         const par::Idx<1> ib{boundaryId};

         const FluxRes fr = characteristicBoundary( surface( mesh.nodes[ip] ),
                                                    dq[ic],q[ic],qb[ib] );
         res[ic]+=fr;
     }
      if( boundaryId==1 )
     {
         const par::Idx<1> ip{mesh.nodes.shape(0)-1};
         const par::Idx<1> ic{mesh.cells.shape(0)-1};
         const par::Idx<1> ib{boundaryId};

         const FluxRes fr = characteristicBoundary( surface( mesh.nodes[ip] ),
                                                    dq[ic],q[ic],qb[ib] );
         res[ic]-=fr;
     }
  }

// =========================================================================================================

/*
 * 1D Higher order flux evaluation with general form and ghost cells at boundaries
 */

   template<LawType Law, ImplementedVarSet SolVarT, ImplementedVarDelta SolDelT, typename HighOrderFlux, floating_point Real>
      requires ConsistentTypes<Law,1,Real,SolVarT,SolDelT>
   void residualCalc( const Mesh<1,Real>&                        mesh,
                      const HighOrderFlux                      hoflux,
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
 * 1D Higher order flux evaluation with general form and periodic boundaries
 */

   template<LawType Law, ImplementedVarSet SolVarT, ImplementedVarDelta SolDelT, typename HighOrderFlux, floating_point Real>
      requires ConsistentTypes<Law,1,Real,SolVarT,SolDelT>
   void residualCalc( const Mesh<1,Real>&                        mesh,
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

