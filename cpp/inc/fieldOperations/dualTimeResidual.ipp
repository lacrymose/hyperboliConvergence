
   template<typename SolutionType>
   void dualTimeResidual( const ODE::Implicit::MultiStep&   bdf,
                          const IdealGas2D::Species&        gas,
                          const Controls::GridControls1D&  grid,
                          const Types::Real                  dt,
                          const std::vector<Array::Array1D<IdealGas2D::VariableSet<SolutionType>>>&   q,
                          const std::vector<Array::Array1D<IdealGas2D::ConservedDelta>>&              r,
                                Array::Array1D<IdealGas2D::ConservedDelta>&                          rm )
  {
      size_t i;
      int j;
   // sanity checks
      assert( rm.size()==grid.n     );
      assert( q.size() ==(size_t)bdf.nsteps );
      assert( r.size() ==(size_t)bdf.nresid );
      for( i=0; i<q.size(); i++ ){ assert( q[i].size()==grid.n ); }
      for( i=0; i<r.size(); i++ ){ assert( r[i].size()==grid.n ); }

      IdealGas2D::ConservedVariables qc;

      rm=0;
   // calculate dualtime source term
      for( i=0; i<grid.n; i++ )
     {
         for( j=0; j<bdf.nsteps; j++ )
        {
            qc = IdealGas2D::ConservedVariables( gas, q[j][i] );
            rm[i]-= bdf.beta[j]*IdealGas2D::ConservedDelta( qc );
        }
     }
      rm/=dt;

   // accumulate flux residuals
      for( j=0; j<bdf.nresid; j++ )
     {
         rm+= bdf.gamma[j]*r[j];
     }

      return;
  }
