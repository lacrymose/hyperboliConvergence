
# pragma once

# include <spatial/boundary/boundaryUpdate.h>
# include <spatial/gradientCalc.h>
# include <spatial/residualCalc.h>
# include <spatial/spectralRadius.h>
# include <spatial/eulerForwardUpdate.h>

# include <conservationLaws/base/base.h>
# include <solutionField/solutionField.h>

# include <lsq/lsq.h>

# include <mesh/mesh.h>

# include <controls.h>

# include <ode.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>
# include <parallalg/parallalg.h>

# include <utils/timing.h>

# include <tuple>
# include <array>
# include <vector>
# include <cassert>

# include <unistd.h>

/*
 * integrates dq/dt = rhs forward in time using an explicit runge kutta scheme
 */
   template<par::execution_policy   Policy,
            LawType                    Law,
            int                       nDim,
            floating_point            Real,
            ImplementedVarSet    SolVarSet,
            typename       SecondOrderFlux,
            typename...      BoundaryConds>
   void integrate( const Policy                                policy,
                   const UnsteadyTimeControls<Real>&     timeControls,
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
      using ResidualArray = par::DualArray<FluxRes,nDim>;

      ResidualArray resTotal(q0.interior.shape());
      std::vector<ResidualArray> resStage = par::vec_of_Arrays<FluxRes,nDim>(rungeKutta.nstages,q0.interior.shape());

   // least squares gradient arrays
      using XMetric = lsq::XMetric<nDim,Real>;
      using QMetric = lsq::QMetric<SolVarSet>;

      using XMetArray = par::DualArray<XMetric,nDim>;
      using QMetArray = par::DualArray<QMetric,nDim>;

   // calculate spatial metrics for least squares
      const XMetric dxdx = xmetrics( policy, mesh.cells );
      QMetric dqdx(mesh.cells.shape());

   // timers
      using FunctionTimer = utils::StopWatchTimer<std::chrono::steady_clock,
                                                  std::chrono::nanoseconds,
                                                  std::chrono::milliseconds>;

      FunctionTimer bcupdate_timer( "bcupdate loop time: " );
      FunctionTimer gradient_timer( "gradient loop time: " );
      FunctionTimer residual_timer( "residual loop time: " );
      FunctionTimer specrads_timer( "specrads loop time: " );
      FunctionTimer rkaccums_timer( "rkaccums loop time: " );
      FunctionTimer eulerfwd_timer( "eulerfwd loop time: " );
      FunctionTimer copyswap_timer( "copyswap func time: " );

      utils::LifetimeTimer timer( "main loop time: " );

      Real t=0;
      for( size_t tstep=0; tstep<timeControls.nTimesteps; tstep++ )
     {
         Real dt{};
         for( unsigned int stg=0; stg<rungeKutta.nstages; stg++ )
        {
         // update boundary conditions
            bcupdate_timer.start();
            boundaryUpdate( mesh,
                            boundaryConds,
                            species,
                            q1 );
            bcupdate_timer.pause();

         // calculate differences
            gradient_timer.start();
            qmetrics( policy,
                      mesh.cells,
                      q1.interior,
                      dqdx );
            gradient_timer.pause();

         // accumulate flux residual
            residual_timer.start();
            residualCalc( policy,
                          flux2,
                          boundaryConds,
                          species,
                          mesh,
                          q1,
                          dxdx,
                          dqdx,
                          resStage[stg] );
            residual_timer.pause();

         // calculate maximum stable timestep for this timestep
            specrads_timer.start();
            if( stg==0 ){ dt = timeControls.cfl/spectralRadius( policy, mesh.cells, resStage[stg] ); }
            specrads_timer.pause();

         // accumulate stage residual
            rkaccums_timer.start();
            rungeKuttaAccumulation( policy,
                                    rungeKutta,
                                    stg,
                                    resStage,
                                    resTotal );
            rkaccums_timer.pause();

         // integrate cell residuals forward by dt and average over cell volume
            eulerfwd_timer.start();
            eulerForwardUpdateGlobal( policy,
                                      mesh.cells,
                                      species,
                                      rungeKutta.beta[stg]*dt,
                                      resTotal,
                                      q0.interior,
                                      q2.interior );
            eulerfwd_timer.pause();

            copyswap_timer.start();
            std::swap( q1,q2 );
            copyswap_timer.pause();
        }
         copyswap_timer.start();
         copy( policy, q0, q1 );
         copyswap_timer.pause();
         t+=dt;
     }
//    std::cout << "physical time elapsed: " << t << "\n";
      std::cout << "\n";
  }

/*
 * Accumulate the total residual for the current runge-kutta stage from the vector of stage residuals
 */
   template<par::execution_policy Policy,
            int                     nDim,
            typename             FluxRes,
            floating_point          Real>
   void rungeKuttaAccumulation( const Policy                                       policy,
                                const ODE::Explicit::RungeKutta<Real>&         rungeKutta,
                                const size_t                                          stg,
                                const std::vector<par::DualArray<FluxRes,nDim>>& resStage,
                                      par::DualArray<FluxRes,nDim>&              resTotal )
  {
   // runge kutta residual accumulation for each cell
      auto rkacc = [&]( const par::DualIdx<nDim>&  idx,
                              FluxRes&         frtotal ) -> void
     {
         for( unsigned int k=0; k<=stg; k++ )
        {
            frtotal.flux+=rungeKutta.alpha[stg][k]*resStage[k](idx).flux;
        }
         return;
     };

      par::fill( policy, resTotal, {} );
      par::for_each_idx( policy, rkacc, resTotal );

      return;
  }

