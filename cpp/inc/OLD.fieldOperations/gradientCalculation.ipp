   template<typename SolutionType>
   void gradientCalculation( const Controls::GridControls1D&   grid,
                             const Array::Array1D< IdealGas2D::VariableSet<   SolutionType > >&  q,
                                   Array::Array1D< IdealGas2D::VariableDelta< SolutionType > >& dq )
  {
      size_t n=grid.n;
      size_t i;
      assert( dq.size() == n+2 );
      assert(  q.size() == n   );

   // internal cells
      for( i=1; i<n; i++ )
     {
         dq[i+1] = q[i]-q[i-1];
     }

   // boundary gradients
      if( grid.boundaryCondition == 'p' )
     {
         dq[  0] = q[n-1] - q[n-2];
         dq[  1] = q[  0] - q[n-1];

         dq[n  ] = q[  0] - q[n-1];
         dq[n+1] = q[  1] - q[  0];
     }
      else if( grid.boundaryCondition == 'a' )
     {
         dq[  0] = 0.;
         dq[  1] = 0.;
         dq[n  ] = 0.;
         dq[n+1] = 0.;
     }

      return;
  }
