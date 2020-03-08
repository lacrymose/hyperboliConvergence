
namespace TimeStepping
{
   template<typename SolutionType, typename Flux, typename Limiter>
   void dualTime( const Controls::TimeSteppingControls&  otime,
                  const ODE::Implicit::MultiStep&          bdf,
                  const Controls::TimeSteppingControls&  itime,
                  const ODE::Explicit::RungeKutta&          rk,
                  const Controls::GridControls1D&         grid,
                  const IdealGas2D::Species&               gas,
                  const Flux                              flux,
                  const Limiter                        limiter,
                        Array::Array1D<IdealGas2D::VariableSet<SolutionType>>& q,
                  const bool print )
  {
      assert( grid.n == q.size() );

      std::vector<Array::Array1D<IdealGas2D::VariableSet<SolutionType>>>    qn(bdf.nsteps);  // outer iteration solutions
      std::vector<Array::Array1D<IdealGas2D::ConservedDelta>>                r(bdf.nresid);   // flux residuals
      std::vector<Array::Array1D<IdealGas2D::ConservedDelta>>               rj(rk.nstages);   // flux residuals
                  Array::Array1D<IdealGas2D::VariableDelta<SolutionType>>   dq;               // solution gradient
                  Array::Array1D<IdealGas2D::ConservedDelta>                rm;               // total dualtime residual
                  Array::Array1D<IdealGas2D::ConservedDelta>                s;                // dual time source term

      IdealGas2D::VariableDelta<SolutionType> du;
      IdealGas2D::ConservedDelta        res,res0;

      Types::Real   lmax,dt,dtau;
      int                    i,j;
      size_t                   k;

      bool  innerConverged;

   // allocate array sizes
      dq.resize( grid.n+2 );
      rm.resize( grid.n );
      s.resize(  grid.n );
      for( i=0; i<bdf.nsteps; i++ ){ qn[i].resize( grid.n ); qn[i]=q; }
      for( i=0; i<bdf.nresid; i++ ){ r[i].resize(  grid.n ); }

   // outer iterations
      for( i=0; i<otime.nt; i++ )
     {
      // q_n residual calculation
         gradientCalculation( grid, qn[0], dq );
         fluxResidual( gas, grid.boundaryCondition,
                       flux, limiter,
                       qn[0],dq, r[bdf.nresid-1],lmax );

      // outer timestep
         if( i==0 ){ dt = otime.cfl/lmax; }

      // inner iterations
         innerConverged=false;
         j=0;
         while( !innerConverged )
        {
         // flux residual
            gradientCalculation( grid, q, dq );
            fluxResidual( gas, grid.boundaryCondition,
                          flux, limiter,
                          q,dq, r[0],lmax );

         // inner timestep
            dtau=itime.cfl/lmax;
            dtau/= ( 1. + bdf.beta[0]*dtau/dt );

         // dualtime residual
            dualTimeResidual( bdf, gas, grid,
                              dt,qn,r, rm );


         // euler update
            eulerForwardNonlinearUpdate( gas, dtau, q,qn[0], rm, res );
            if( j==0 ){ res0=res; }
            q=qn[0];

         // print inner residual
            std::cout << std::left << std::setw(4) << i;
            std::cout << std::left << std::setw(4) << j;
            std::cout << res << std::endl;

         // test convergence
            j++;
            innerConverged= true;
            for( k=0; k<4; k++ )
           {
               innerConverged = innerConverged and ( res[k]/res0[k] < itime.convergenceDrop );
           }
            innerConverged = innerConverged or ( j > itime.nt );
        }

      // print outer residual
         residual( gas, qn[0],qn[1], res );

         std::cout << std::endl;
         std::cout << std::left << std::setw(8) << i;
         std::cout << res << std::endl;
         std::cout << std::endl;

      // rearrange
         for( k=bdf.nsteps-1; k>0; k-- ){ qn[k]=qn[k-1]; }
     }
      return;
  }
}

