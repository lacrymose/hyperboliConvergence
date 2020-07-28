
namespace TimeStepping
{
   template<typename SolutionType, typename Flux, typename Limiter>
   void rungeKutta( const Controls::TimeSteppingControls& time,
                    const ODE::Explicit::RungeKutta&        rk,
                    const Controls::GridControls1D&       grid,
                    const IdealGas2D::Species&             gas,
                    const Flux                            flux,
                    const Limiter                      limiter,
                          Array::Array1D< IdealGas2D::VariableSet< SolutionType > >& q,
                    const bool print )
  {
      assert( grid.n == q.size() );

      Array::Array1D<IdealGas2D::VariableDelta<SolutionType>>   dq; // solution gradient
      Array::Array1D<IdealGas2D::VariableSet<  SolutionType>>   q1; // solution at end of each rk stage
      Array::Array1D<IdealGas2D::ConservedDelta>                rt; // total residual for each rk stage

      std::vector<Array::Array1D<IdealGas2D::ConservedDelta>>  r(rk.nstages,0); // residual evaluation at each rk stage

      IdealGas2D::ConservedDelta                             res;
      Types::Real lmax,dt;

      assert( grid.n == q.size() );

      dq.resize( grid.n+2 );
      q1.resize( grid.n );
      rt.resize( grid.n );

      for( int j=0; j<rk.nstages; j++ ){ r[j].resize( grid.n ); }

      q1=q;
      for( int i=0; i<time.nt; i++ )
     {
         lmax=-1;

         for( int j=0; j<rk.nstages; j++ )
        {
         // flux evaluation
            gradientCalculation( grid, q1, dq );
            fluxResidual( gas, grid.boundaryCondition,
                          flux, limiter,
                          q1,dq, r[j],lmax );

            if( j==0 ){ dt = time.cfl/lmax; }

         // accumulate stage residual
            rt=0;
            for( int k=0; k<=j; k++ )
           {
               rt+=rk.alpha[j][k]*r[k];
           }

         // update along stage-scaled timestep
            eulerForwardNonlinearUpdate( gas, rk.beta[j]*dt, q, q1, rt, res );
        }
         q=q1;

         if( print ){ std::cout << res << std::endl; }
     }

      return;
  }
}

/*
 * integrates dq/dt = rhs forward in time using an explicit runge kutta scheme
 */
   template<int                   nDim,
           ImplementedVarSet   SolVarT,
           typename    ResidualFunctor,
           typename    TimestepFunctor,
           typename      UpdateFunctor>
   void integrate( const UnsteadyTimeControls&        utc,
                   const ODE::ExplicitRungeKutta&      rk,
                   const Mesh<nDim,Real>&            mesh,
                         par::Array<SolVarSet,nDim>&   q0,
                         RHSFunctor&                  rhs,
                         TimestepFunctor&            delt,
                         UpdateFunctor&            update )
  {
   // check sizes match
      assert( q0.shape() == mesh.shape() );

      const par::Shape<nDim> shape = q0.shape();

      par::Array<FluxRes,nDim>               resTotal(shape);
      std::vector<par::Array<FluxRes,nDim>>  resStage = par::vec_of_Arrays<FluxRes,nDim>(rk.nstages,shape);

      Real t=0;
      for( size_t tstep=0; tstep<tc.nt; tstep++ )
     {
         Real dt;
         for( int stg=0; stg<rk.nstages; stg++ )
        {
            resStage[stg] = rhs( mesh, t, q1 );

            if( stg==0 ){ dt = delt( mesh, resStage[0] ); }

            resTotal = weightedSum( stg+1, rk.alpha[stg], resStage );

            q2 = UpdateFunctor( mesh, rk.beta[stg]*dt, resTotal, q0 );
            std::swap( q1, q2 );
        }
         par::copy( q0, q1 );
         t+=dt;
     }
  }

