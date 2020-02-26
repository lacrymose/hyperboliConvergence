
namespace TimeStepping
{
   template<typename Flux, typename SolutionType>
   void ExplicitEuler( const Controls::TimeSteppingControls& time,
                       const Controls::GridControls1D&       grid,
                       const IdealGas2D::Species&             gas,
                       const Flux                            flux,
                             Array::Array1D< IdealGas2D::VariableSet< SolutionType > >& q,
                       const bool print )
  {

      Array::Array1D<IdealGas2D::VariableSet<SolutionType>>   q1;
      Array::Array1D<IdealGas2D::ConservedDelta>               r;
      IdealGas2D::ConservedDelta                             res;
      Types::Real lmax,dt;

      assert( grid.n == q.size() );

      q1.resize( grid.n );
      r.resize(  grid.n );

      for( int i=0; i<time.nt; i++ )
     {
         lmax=-1;

         fluxResidual( gas, grid.boundaryCondition, flux, q,r, lmax );

         dt = time.cfl/lmax;

         eulerForwardNonlinearUpdate( gas, dt, q, q1, r, res );
         q=q1;

         if( print ){ std::cout << res << std::endl; }
     }

      return;
  }
}
