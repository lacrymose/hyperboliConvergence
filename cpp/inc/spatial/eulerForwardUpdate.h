
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <vector>

# include <cassert>

// overload with return value (must be used to construct vector to use RVO)
   template<LawType Law, int nDim, ImplementedVarSet SolVarT, floating_point Real>
      requires ConsistentTypes<Law,nDim,Real,SolVarT>
   par::Array<SolVarT,nDim> eulerForwardUpdate( const par::Array<geom::Volume<nDim,Real>,nDim>& cells,
                                                const Species<Law,Real>&                      species,
                                                const Real                                        cfl,
                                                const Real                                       lmax,
                                                const par::Array<FluxResult<Law,nDim,Real>,nDim>&   r,
                                                const par::Array<SolVarT,nDim>&                    q0 )
  {
      par::Array<SolVarT,nDim> q1(q0.shape());
      eulerForwardUpdate( cells, species, cfl, lmax, r, q0, q1 );
      return q1;
  }

   template<LawType Law, int nDim, ImplementedVarSet SolVarT, floating_point Real>
      requires ConsistentTypes<Law,nDim,Real,SolVarT>
   void eulerForwardUpdate( const par::Array<geom::Volume<nDim,Real>,nDim>& cells,
                            const Species<Law,Real>&                      species,
                            const Real                                        cfl,
                            const Real                                       lmax,
                            const par::Array<FluxResult<Law,nDim,Real>,nDim>&   r,
                            const par::Array<SolVarT,nDim>&                    q0,
                                  par::Array<SolVarT,nDim>&                    q1 )
  {
   // conserved variables/deltas needed for correct shock speeds
      constexpr BasisType<Law> ConservedBasis = BasisType<Law>::Conserved;
      using ConsVarT = VariableSet<  Law,nDim,ConservedBasis,Real>;
      using ConsDelT = VariableDelta<Law,nDim,ConservedBasis,Real>;

   // check mesh sizes match
      assert( cells.shape() == q1.shape() );
      assert( cells.shape() == q0.shape() );
      assert( cells.shape() ==  r.shape() );

//    const size_t nc = cells.dims[0];

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
      // solv -> consv + increment -> solv
      const auto update = [&consvar, &solvar, dt]
                          ( const SolVarT& v0,
                            const geom::Volume<nDim,Real>&   cell,
                            const FluxResult<Law,nDim,Real>&  res ) -> SolVarT
                          {
                             const ConsDelT dvc = res.flux*(dt/cell.volume);
                             return solvar( consvar( v0 ) + dvc );
                          };

      par::transform( update,
                      q1,
                      q0, cells, r );
      return;
  }
