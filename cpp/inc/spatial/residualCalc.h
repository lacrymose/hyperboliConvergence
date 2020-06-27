
# pragma once

# include <limiters/limiter.h>

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <mdarray/mdarray.h>

# include <vector>

# include <cassert>

// overload with return value
   template<LawType Law, ImplementedVarSet SolVarT, FluxFunctor<Law> Flux, floating_point Real>
      requires ConsistentTypes<Law,1,Real,SolVarT>
   std::vector<FluxResult<Law,1,Real>> residualCalc( const std::vector<geom::Volume<1,Real>>& cells,
                                                     const std::vector<geom::Point< 1,Real>>& nodes,
                                                     const Flux flux,
                                                     const Species<Law,Real>& species,
                                                     const SolVarT& qlb,
                                                     const std::vector<SolVarT>& q )
  {
      using FluxRes = FluxResult<Law,1,Real>;
      std::vector<FluxRes> res(q.size());
      residualCalc( cells, nodes, flux, species, qlb, q, res );
      return res;
  }

/*
 * First order flux evaluation
 */
   template<LawType Law, ImplementedVarSet SolVarT, FluxFunctor<Law> Flux, floating_point Real>
      requires ConsistentTypes<Law,1,Real,SolVarT>
   void residualCalc( const std::vector<geom::Volume<1,Real>>& cells,
                      const std::vector<geom::Point< 1,Real>>& nodes,
                      const Flux flux,
                      const Species<Law,Real>& species,
                      const SolVarT& qlb,
                      const std::vector<SolVarT>& q,
                            std::vector<FluxResult<Law,1,Real>>& res )
  {
      using FluxRes = FluxResult<Law,1,Real>;
      using Surface = geom::Surface<1,Real>;

   // check mesh sizes match
      const size_t nc = cells.size();

      assert( nc+1 == nodes.size() );
      assert( nc   == res.size() );
      assert( nc   ==   q.size() );

   // reset residual
      for( FluxRes& fr : res ){ fr = FluxRes{}; }

   // accumulate cell residual contributions from the flux across each face
      for( size_t i=0; i<nc-1; i++ )
     {
         const SolVarT ql=q[i];
         const SolVarT qr=q[i+1];

         const Surface face = surface( nodes[i] );

         const FluxRes fr = flux( species, face,  ql, qr );

         res[i]  -=fr;
         res[i+1]+=fr;
     }

   // left boundary : dirichlet condition at initial state
     {
         const size_t i=0;
         const SolVarT qr=q[i];
         const Surface face = surface( nodes[i] );
         const FluxRes fr = flux( species, face, qlb, qr );
         res[i]+=fr;
     }

   // right boundary : outflow, upwind flux
     {
         const size_t i=nc-1;
         const SolVarT ql=q[i];
         const Surface face = surface( nodes[i] );
         const FluxRes fr = flux( species, face, ql, ql );
         res[i]-=fr;
     }
      return;
  }

// overload with return value
   template<LawType Law, ImplementedVarSet SolVarT, ImplementedVarDelta SolDelT, typename Limiter, FluxFunctor<Law> Flux, floating_point Real>
      requires ConsistentTypes<Law,1,Real,SolVarT,SolDelT>
   std::vector<FluxResult<Law,1,Real>> residualCalc( const std::vector<geom::Volume<1,Real>>& cells,
                                                     const std::vector<geom::Point< 1,Real>>& nodes,
                                                     const Flux                    flux,
                                                     const Limiter              limiter,
                                                     const Species<Law,Real>&   species,
                                                     const SolVarT&                 qlb,
                                                     const std::vector<SolVarT>&      q,
                                                     const std::vector<SolDelT>&     dq )
  {
      using FluxRes = FluxResult<Law,1,Real>;
      std::vector<FluxRes> res(q.size());
      residualCalc( cells, nodes, flux, limiter, species, qlb, q,dq, res );
      return res;
  }

/*
 * Higher order flux evaluation
 */
   template<LawType Law, ImplementedVarSet SolVarT, ImplementedVarDelta SolDelT, typename Limiter, FluxFunctor<Law> Flux, floating_point Real>
      requires ConsistentTypes<Law,1,Real,SolVarT,SolDelT>
   void residualCalc( const std::vector<geom::Volume<1,Real>>& cells,
                      const std::vector<geom::Point< 1,Real>>& nodes,
                      const Flux                    flux,
                      const Limiter              limiter,
                      const Species<Law,Real>&   species,
                      const SolVarT&                 qlb,
                      const std::vector<SolVarT>&      q,
                      const std::vector<SolDelT>&     dq,
                            std::vector<FluxResult<Law,1,Real>>& res )
  {
      using FluxRes = FluxResult<Law,1,Real>;
      using Surface = geom::Surface<1,Real>;

   // check mesh sizes match
      const size_t nc = cells.size();

      assert( nc+1 == nodes.size() );
      assert( nc   == res.size() );
      assert( nc   ==  dq.size() );
      assert( nc   ==   q.size() );

   // limit all elements in vardelta
      const auto limvar = [&limiter]( const SolDelT& dq0, const SolDelT& dq1 ) -> SolDelT
     {
         SolDelT dqlim;
         for( int i=0; i<SolDelT::N; i++ )
        {
            dqlim[i] = limiter( dq0[i],dq1[i] );
        }
         return dqlim;
     };

   // reset residual
      for( FluxRes& fr : res ){ fr = FluxRes{}; }

   // accumulate cell residual contributions from the flux across each face
      for( size_t i=0; i<nc-1; i++ )
     {
      // left/right cell values
         const SolVarT ql0=q[i];
         const SolVarT qr0=q[i+1];

      // central, and left/right bias differences
         const SolDelT dqc = qr0 - ql0;
         const SolDelT dql = dq[i  ] - dqc;
         const SolDelT dqr = dq[i+1] - dqc;

         const SolDelT slopel = limvar( dqc, dql );
         const SolDelT sloper = limvar( dqc, dqr );

      // interface values
         const SolVarT ql1 = ql0 + 0.5*slopel;
         const SolVarT qr1 = qr0 - 0.5*sloper;

         const Surface face = surface( nodes[i] );

         const FluxRes fr = flux( species, face,  ql1, qr1 );

         res[i]  -=fr;
         res[i+1]+=fr;
     }

   // first order boundaries
   // left boundary : dirichlet condition at initial state
     {
         const size_t i=0;
         const SolVarT qr=q[i];
         const Surface face = surface( nodes[i] );
         const FluxRes fr = flux( species, face, qlb, qr );
         res[i]+=fr;
     }

   // right boundary : outflow, upwind flux
     {
         const size_t i=nc-1;
         const SolVarT ql=q[i];
         const Surface face = surface( nodes[i] );
         const FluxRes fr = flux( species, face, ql, ql );
         res[i]-=fr;
     }

      return;
  }

/*
 * Higher order flux evaluation with general form
 */

   template<LawType Law, ImplementedVarSet SolVarT, ImplementedVarDelta SolDelT, typename HighOrderFlux, floating_point Real>
      requires ConsistentTypes<Law,1,Real,SolVarT,SolDelT>
   void residualCalc( const MDArray<geom::Volume<1,Real>,1>& cells,
                      const MDArray<geom::Point< 1,Real>,1>& nodes,
                      const HighOrderFlux                   hoflux,
                      const SolVarT&                           qlb,
                      const MDArray<SolVarT,1>&                  q,
                      const MDArray<SolDelT,1>&                 dq,
                            MDArray<FluxResult<Law,1,Real>,1>& res )
  {
      using FluxRes = FluxResult<Law,1,Real>;
      using Surface = geom::Surface<1,Real>;

   // check mesh sizes match

      assert( cells.dims == res.dims );
      assert( cells.dims ==  dq.dims );
      assert( cells.dims ==   q.dims );
      assert( nodes.dims == Dims<1>{cells.dims[0]+1} );

      const size_t nc = cells.dims[0];

   // reset residual
      for( FluxRes& fr : res.elems ){ fr = FluxRes{}; }

   // accumulate cell residual contributions from the flux across each face
      for( size_t i=0; i<nc-1; i++ )
     {
      // left/right cell values
         const SolVarT ql=q[{i}];
         const SolVarT qr=q[{i+1}];

      // left,right and central bias differences
         const SolDelT dqc = qr - ql;
         const SolDelT dql = dq[{i  }] - dqc;
         const SolDelT dqr = dq[{i+1}] - dqc;

         const Surface face = surface( nodes[{i}] );

         const FluxRes fr = hoflux( face, dql,dqc,dqr, ql,qr );

         res[{i}]  -=fr;
         res[{i+1}]+=fr;
     }

   // first order boundaries
   // left boundary : dirichlet condition at initial state
     {
         const size_t i=0;
         const SolVarT ql=qlb;
         const SolVarT qr=q[{i}];

      // zero gradients
         const SolDelT dqc = SolDelT{};
         const SolDelT dql = SolDelT{};
         const SolDelT dqr = SolDelT{};

         const Surface face = surface( nodes[{i}] );
         const FluxRes fr = hoflux( face, dql,dqc,dqr, ql,qr );
         res[{i}]+=fr;
     }

   // right boundary : outflow, upwind flux
     {
         const size_t i=nc-1;
         const SolVarT ql=q[{i}];
         const SolVarT qr=q[{i}];

      // zero gradients
         const SolDelT dqc = SolDelT{};
         const SolDelT dql = SolDelT{};
         const SolDelT dqr = SolDelT{};

         const Surface face = surface( nodes[{i}] );
         const FluxRes fr = hoflux( face, dql,dqc,dqr, ql,qr );
         res[{i}]-=fr;
     }

      return;
  }

