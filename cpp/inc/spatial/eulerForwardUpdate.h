
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <vector>

# include <cassert>

// overload with return value (must be used to construct vector to use RVO)
   template<LawType Law, int nDim, ImplementedVarSet SolVarT, floating_point Real>
      requires ConsistentTypes<Law,nDim,Real,SolVarT>
   std::vector<SolVarT> eulerForwardUpdate( const std::vector<geom::Volume<nDim,Real>>& cells,
                                            const Species<Law,Real>& species,
                                            const Real cfl,
                                            const Real lmax,
                                            const std::vector<FluxResult<Law,nDim,Real>>& r,
                                            const std::vector<SolVarT>& q0 )
  {
      std::vector<SolVarT> q1(q0.size());
      eulerForwardUpdate( cells, species, cfl, lmax, r, q0, q1 );
      return q1;
  }

   template<LawType Law, int nDim, ImplementedVarSet SolVarT, floating_point Real>
      requires ConsistentTypes<Law,nDim,Real,SolVarT>
   void eulerForwardUpdate( const std::vector<geom::Volume<nDim,Real>>& cells,
                            const Species<Law,Real>& species,
                            const Real cfl,
                            const Real lmax,
                            const std::vector<FluxResult<Law,nDim,Real>>& r,
                            const std::vector<SolVarT>& q0,
                                  std::vector<SolVarT>& q1 )
  {
      constexpr BasisType<Law> ConservedBasis = BasisType<Law>::Conserved;
      using ConsVarT = VariableSet<  Law,nDim,ConservedBasis,Real>;
      using ConsDelT = VariableDelta<Law,nDim,ConservedBasis,Real>;

   // check mesh sizes match
      const size_t nc = cells.size();

      assert( nc == q1.size() );
      assert( nc == q0.size() );
      assert( nc ==  r.size() );

   // timestep
      const Real dt = cfl/lmax;

      for( size_t i=0; i<nc; i++ )
     {
         const Real d = dt/cells[i].volume;

         const ConsDelT dqc = d*r[i].flux;

      // transform to conserved variables, increment, and transform back
         q1[i] = set2Set<SolVarT>( species,
                                   set2Set<ConsVarT>( species,
                                                      q0[i] ) + dqc );
     }
      return;
  }
