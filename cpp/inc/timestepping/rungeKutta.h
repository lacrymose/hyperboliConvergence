
# pragma once

# include <spatial/boundary/boundaryUpdate.h>
# include <spatial/gradientCalc.h>
# include <spatial/residualCalc.h>
# include <spatial/spectralRadius.h>
# include <spatial/eulerForwardUpdate.h>

# include <conservationLaws/base/base.h>
# include <solutionField/solutionField.h>

# include <mesh/mesh.h>

# include <controls.h>

# include <ode.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <tuple>
# include <array>
# include <vector>
# include <cassert>

/*
 * integrates dq/dt = rhs forward in time using an explicit runge kutta scheme
 */
   template<LawType                    Law,
            int                       nDim,
            floating_point            Real,
            ImplementedVarSet    SolVarSet,
            typename       SecondOrderFlux,
            typename...      BoundaryConds>
   void integrate( const UnsteadyTimeControls<Real>&     timeControls,
                   const ODE::Explicit::RungeKutta<Real>&  rungeKutta,
                   const SecondOrderFlux&                       flux2,
                   const std::tuple<BoundaryConds...>   boundaryConds,
                   const Species<Law,Real>&                   species,
                   const Mesh<nDim,Real>&                        mesh,
                         SolutionField<SolVarSet,nDim>&            q0 )
  {
   // check sizes match
      assert( q0.interior.shape() == mesh.cells.shape() );

   // spare solution arrays for rk/timestepping iterations
      // q0 is solution at beginning of current timestep
      // q1 is solution at beginning of current rk stage
      // q2 is working array to save update into at end of rk stage
      SolutionField<SolVarSet,nDim> q1 = copy(q0);
      SolutionField<SolVarSet,nDim> q2 = copy(q0);

   // residual arrays
      using FluxRes = FluxResult<Law,nDim,Real>;

      par::Array<FluxRes,nDim>  resTotal(q0.interior.shape());
      std::vector<par::Array<FluxRes,nDim>>  resStage = par::vec_of_Arrays<FluxRes,nDim>(rungeKutta.nstages,q0.interior.shape());

   // gradient arrays
      using SolVarDel = vardelta_t<SolVarSet>;
      using SolVarGrad = std::array<SolVarDel,nDim>;

      par::Array<SolVarGrad,nDim> dq( q0.interior.shape(), SolVarGrad{} );

      for( size_t tstep=0; tstep<timeControls.nTimesteps; tstep++ )
     {
         Real lmax{};
         for( unsigned int stg=0; stg<rungeKutta.nstages; stg++ )
        {
         // update boundary conditions
            boundaryUpdate( mesh,
                            boundaryConds,
                            species,
                            q1 );

         // calculate differences
            gradientCalc( mesh,
                          q1,
                          dq );

         // accumulate flux residual
            residualCalc( mesh,
                          flux2,
                          boundaryConds,
                          species,
                          q1,
                          dq,
                          resStage[stg] );

         // calculate maximum stable timestep for this timestep
            if( stg==0 ){ lmax = spectralRadius( mesh.cells, resStage[stg] ); }

         // accumulate stage residual
            par::fill( resTotal, {} );

         // runge-kutta accumulation closure
            rungeKuttaAccumulation( rungeKutta,
                                    stg,
                                    resStage,
                                    resTotal );

         // integrate cell residuals forward by dt and average over cell volume
            eulerForwardUpdate( mesh.cells,
                                species,
                                rungeKutta.beta[stg]*timeControls.cfl/lmax,
                                resTotal,
                                q0.interior,
                                q2.interior );
            std::swap( q1,q2 );
        }
         copy( q0, q1 );
     }
  }

/*
 * Accumulate the total residual for the current runge-kutta stage from the vector of stage residuals
 */
   template<int              nDim,
            typename      FluxRes,
            floating_point   Real>
   void rungeKuttaAccumulation( const ODE::Explicit::RungeKutta<Real>&    rungeKutta,
                                const size_t                                     stg,
                                const std::vector<par::Array<FluxRes,nDim>>& resStage,
                                      par::Array<FluxRes,nDim>&              resTotal )
  {
   // runge kutta residual accumulation for each cell
      auto rkacc = [&]( const par::Idx<nDim>&  idx,
                              FluxRes&     frtotal ) -> void
     {
         for( unsigned int k=0; k<=stg; k++ )
        {
            frtotal.flux+=rungeKutta.alpha[stg][k]*resStage[k][idx].flux;
        }
         return;
     };
      par::for_each_idx( rkacc, resTotal );

      return;
  }

