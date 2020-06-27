
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <vector>

# include <cassert>

// overload with return value (must be used to construct vector to use RVO)
   template<LawType Law, int nDim, ImplementedVarSet SolVarT, floating_point Real>
      requires ConsistentTypes<Law,nDim,Real,SolVarT>
   MDArray<SolVarT,1> eulerForwardUpdate( const MDArray<geom::Volume<nDim,Real>,1>& cells,
                                          const Species<Law,Real>& species,
                                          const Real cfl,
                                          const Real lmax,
                                          const MDArray<FluxResult<Law,nDim,Real>,1>& r,
                                          const MDArray<SolVarT,1>& q0 )
  {
      MDArray<SolVarT,1> q1(q0.dims);
      eulerForwardUpdate( cells, species, cfl, lmax, r, q0, q1 );
      return q1;
  }

   template<LawType Law, int nDim, ImplementedVarSet SolVarT, floating_point Real>
      requires ConsistentTypes<Law,nDim,Real,SolVarT>
   void eulerForwardUpdate( const MDArray<geom::Volume<nDim,Real>,1>& cells,
                            const Species<Law,Real>& species,
                            const Real cfl,
                            const Real lmax,
                            const MDArray<FluxResult<Law,nDim,Real>,1>& r,
                            const MDArray<SolVarT,1>& q0,
                                  MDArray<SolVarT,1>& q1 )
  {
   // conserved variables/deltas needed for correct shock speeds
      constexpr BasisType<Law> ConservedBasis = BasisType<Law>::Conserved;
      using ConsVarT = VariableSet<  Law,nDim,ConservedBasis,Real>;
      using ConsDelT = VariableDelta<Law,nDim,ConservedBasis,Real>;

   // check mesh sizes match
      assert( cells.dims == q1.dims );
      assert( cells.dims == q0.dims );
      assert( cells.dims ==  r.dims );

      const size_t nc = cells.dims[0];

   // timestep
      const Real dt = cfl/lmax;

   // transformation to conservative variables
      const auto consvar = [&species]( const SolVarT& qs ) -> ConsVarT
     {
         return set2Set<ConsVarT>( species, qs );
     };

   // transformation to solution variables
      const auto solvar = [&species]( const ConsVarT& qc ) -> SolVarT
     {
         return set2Set<SolVarT>( species, qc );
     };

   // new = old + dt*residual/vol
      for( size_t i=0; i<nc; i++ )
     {
         const Real d = dt/cells[{i}].volume;

         const ConsDelT dqc = d*r[{i}].flux;

      // solv -> consv + increment -> solv
         q1[{i}] = solvar( consvar( q0[{i}] ) + dqc );
     }
      return;
  }
