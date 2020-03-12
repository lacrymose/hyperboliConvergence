
   template<typename T>
   void centralExplicitSmoothing( const Controls::GridControls1D&  grid,
                                  const Types::Real               alpha,
                                  const Array::Array1D<T>&           r0,
                                        Array::Array1D<T>&           r1 )
  {
      size_t i,n=grid.n;
      assert( r0.size() == n );
      assert( r1.size() == n );

      Types::Real beta,gamma;
      auto  dr=r0[0]-r0[1];

      beta = 0.25*( alpha*alpha - 1. );

      gamma = 0.25*( (1.+4.*beta)/alpha - 1. );

      r1=r0;

      for( i=1; i<n; i++ )
     {
         dr = gamma*( r0[i] - r0[i-1] );
         r1[i-1] -= dr;
         r1[i  ] += dr;
     }

      if( grid.boundaryCondition=='p' )
     {
         dr = gamma*( r0[0] - r0[n-1] );
         r1[n-1] -= dr;
         r1[0  ] += dr;
     }

      return;
  }

   template<typename T>
   void centralExplicitSmoothing2( const Controls::GridControls1D&  grid,
                                   const Types::Real               alpha,
                                   const Array::Array1D<T>&           r0,
                                         Array::Array1D<T>&           r1 )
  {
      size_t i,n=grid.n;
      assert( r0.size() == n );
      assert( r1.size() == n );

      Types::Real beta,gamma;
      T  dr;

      beta = 0.25*( alpha*alpha - 1. );

      gamma = 0.25*( (1.+4.*beta)/alpha - 1. );

      r1=r0;

      for( i=1; i<n-1; i++ )
     {
         r1[i] = r0[i] + gamma*( r0[i-1] + r0[i+1] );
         r1[i]/= 1.+2.*gamma;
     }

      if( grid.boundaryCondition=='p' )
     {
         r1[  0] = r0[  0] + gamma*( r0[n-1] + r0[  1] );
         r1[n-1] = r0[n-1] + gamma*( r0[n-2] + r0[  0] );

         r1[  0]/= 1.+2.*gamma;
         r1[n-1]/= 1.+2.*gamma;
     }
      else if( grid.boundaryCondition=='a' )
     {
         r1[  0] = r0[  0] + gamma*(           r0[  1] );
         r1[n-1] = r0[n-1] + gamma*( r0[n-2]           );

         r1[  0]/= 1.+1.*gamma;
         r1[n-1]/= 1.+1.*gamma;
     }

      return;
  }

