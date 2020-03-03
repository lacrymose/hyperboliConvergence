
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
      Array::Array1D<IdealGas2D::VariableSet<SolutionType>>     q1; // solution at end of each rk stage
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
                          flux,limiter,
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
