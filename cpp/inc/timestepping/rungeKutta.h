
# pragma once

# include <spatial/spectralRadius.h>
# include <spatial/rungeKuttaAccumulation.h>
# include <spatial/eulerForwardUpdate.h>

# include <conservationLaws/base/base.h>

# include <mesh/mesh.h>

# include <controls.h>

# include <ode.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <vector>
# include <cassert>

/*
 * integrates dq/dt = rhs forward in time using an explicit runge kutta scheme
 */
   template<LawType                 Law,
            int                    nDim,
            floating_point         Real,
            ImplementedVarSet SolVarSet,
            typename    ResidualFunctor>
   void integrate( const UnsteadyTimeControls<Real>&     timeControls,
                   const ODE::Explicit::RungeKutta<Real>&  rungeKutta,
                   const Species<Law,Real>&                   species,
                   const Mesh<nDim,Real>&                        mesh,
                         par::Array<SolVarSet,nDim>&               q0,
                         ResidualFunctor&                         rhs )
  {
   // check sizes match
      assert( q0.shape() == mesh.shape );

      par::Array<SolVarSet,nDim> q1 = par::copy(q0);

      using FluxRes = FluxResult<Law,nDim,Real>;

      par::Array<FluxRes,nDim>  resTotal(q0.shape());
      std::vector<par::Array<FluxRes,nDim>>  resStage = par::vec_of_Arrays<FluxRes,nDim>(rungeKutta.nstages,q0.shape());

      Real t=0;
      for( size_t tstep=0; tstep<timeControls.nTimesteps; tstep++ )
     {
         Real dt{};
         for( unsigned int stg=0; stg<rungeKutta.nstages; stg++ )
        {
            resStage[stg] = rhs( mesh, t, q1 );

            if( stg==0 ){ dt = timeControls.cfl/spectralRadius( mesh.cells, resStage[0] ); }

            resTotal = rungeKuttaAccumulation( rungeKutta, stg, resStage );

            const Real dt_stage = rungeKutta.beta[stg]*dt;
            q1 = eulerForwardUpdate( mesh.cells, species, dt_stage, resTotal, q0 );
        }
         q0 = par::copy( q1 );
         t+=dt;
     }
  }

