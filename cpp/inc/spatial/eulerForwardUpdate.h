
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>
# include <parallalg/parallalg.h>

# include <vector>

# include <cassert>

// overload with return value (must be used to construct vector to use RVO)
   template<par::execution_policy Policy,
            LawType                  Law,
            int                     nDim,
            ImplementedVarSet    SolVarT,
            floating_point          Real>
      requires ConsistentTypes<Law,nDim,Real,SolVarT>
   par::Array<SolVarT,nDim> eulerForwardUpdateGlobal( const Policy                                       policy,
                                                      const par::DualArray<geom::Volume<nDim,Real>,nDim>& cells,
                                                      const Species<Law,Real>&                          species,
                                                      const Real                                             dt,
                                                      const par::DualArray<FluxResult<Law,nDim,Real>,nDim>&   r,
                                                      const par::DualArray<SolVarT,nDim>&                    q0 )
  {
      par::DualArray<SolVarT,nDim> q1(q0.shape());
      eulerForwardUpdateGlobal( policy, cells, species, dt, r, q0, q1 );
      return q1;
  }

/*
 * One step of euler forward using residual vector r and a global timestep dt
 */
   template<par::execution_policy Policy,
            LawType                  Law,
            int                     nDim,
            ImplementedVarSet    SolVarT,
            floating_point          Real>
      requires ConsistentTypes<Law,nDim,Real,SolVarT>
   void eulerForwardUpdateGlobal( const Policy                                       policy,
                                  const par::DualArray<geom::Volume<nDim,Real>,nDim>& cells,
                                  const Species<Law,Real>&                          species,
                                  const Real                                             dt,
                                  const par::DualArray<FluxResult<Law,nDim,Real>,nDim>&   r,
                                  const par::DualArray<SolVarT,nDim>&                    q0,
                                        par::DualArray<SolVarT,nDim>&                    q1 )
  {
   // conserved variables/deltas needed for correct shock speeds
      constexpr BasisType<Law> ConservedBasis = BasisType<Law>::Conserved;
      using ConsVarT = VariableSet<  Law,nDim,ConservedBasis,Real>;
      using ConsDelT = VariableDelta<Law,nDim,ConservedBasis,Real>;

   // check mesh sizes match
      assert( cells.shape() == q1.shape() );
      assert( cells.shape() == q0.shape() );
      assert( cells.shape() ==  r.shape() );

   // transformation to conservative variables
      const auto consvar = [&species]
                           ( const SolVarT& qs ) -> ConsVarT
     {
         return set2Set<ConsVarT>( species, qs );
     };

   // transformation to solution variables
      const auto solvar = [&species]
                          ( const ConsVarT& qc ) -> SolVarT
     {
         return set2Set<SolVarT>( species, qc );
     };

   // new = old + dt*residual/vol
      // solv -> (consv + increment) -> solv
      const auto update = [&consvar, &solvar, dt]
                          ( const SolVarT&                     v0,
                            const geom::Volume<nDim,Real>&   cell,
                            const FluxResult<Law,nDim,Real>&  res ) -> SolVarT
     {
         const ConsDelT dvc = res.flux*(dt/cell.volume);
         return solvar( consvar( v0 ) + dvc );
     };

      par::transform( policy,
                      update,
                      q1,
                      q0, cells, r );
      return;
  }



// overload with return value (must be used to construct vector to use RVO)
   template<par::execution_policy Policy,
            LawType                  Law,
            int                     nDim,
            ImplementedVarSet    SolVarT,
            floating_point          Real>
      requires ConsistentTypes<Law,nDim,Real,SolVarT>
   par::Array<SolVarT,nDim> eulerForwardUpdateLocal( const Policy                                       policy,
                                                     const Species<Law,Real>&                          species,
                                                     const Real                                            cfl,
                                                     const par::DualArray<FluxResult<Law,nDim,Real>,nDim>&   r,
                                                     const par::DualArray<SolVarT,nDim>&                    q0 )
  {
      par::DualArray<SolVarT,nDim> q1(q0.shape());
      eulerForwardUpdateLocal( policy, species, cfl, r, q0, q1 );
      return q1;
  }

/*
 * One step of euler forward using residual vector r and local timestepping timestep with cfl
 */
   template<par::execution_policy Policy,
            LawType                  Law,
            int                     nDim,
            ImplementedVarSet    SolVarT,
            floating_point          Real>
      requires ConsistentTypes<Law,nDim,Real,SolVarT>
   void eulerForwardUpdateLocal( const Policy                                       policy,
                                 const Species<Law,Real>&                          species,
                                 const Real                                            cfl,
                                 const par::DualArray<FluxResult<Law,nDim,Real>,nDim>&   r,
                                 const par::DualArray<SolVarT,nDim>&                    q0,
                                       par::DualArray<SolVarT,nDim>&                    q1 )
  {
   // conserved variables/deltas needed for correct shock speeds
      constexpr BasisType<Law> ConservedBasis = BasisType<Law>::Conserved;
      using ConsVarT = VariableSet<  Law,nDim,ConservedBasis,Real>;
      using ConsDelT = VariableDelta<Law,nDim,ConservedBasis,Real>;

   // check mesh sizes match
      assert( r.shape() == q1.shape() );
      assert( r.shape() == q0.shape() );

   // transformation to conservative variables
      const auto consvar = [&species]
                           ( const SolVarT& qs ) -> ConsVarT
     {
         return set2Set<ConsVarT>( species, qs );
     };

   // transformation to solution variables
      const auto solvar = [&species]
                          ( const ConsVarT& qc ) -> SolVarT
     {
         return set2Set<SolVarT>( species, qc );
     };

   // new = old + dt*residual/vol
      // solv -> (consv + increment) -> solv
      const auto update = [consvar, solvar, cfl]
                          ( const SolVarT&                     v0,
                            const FluxResult<Law,nDim,Real>&  res ) -> SolVarT
     {
         const ConsDelT dvc = res.flux*(cfl/res.lambda);
         return solvar( consvar( v0 ) + dvc );
     };

      par::transform( policy,
                      update,
                      q1,
                      q0, r );
      return;
  }
