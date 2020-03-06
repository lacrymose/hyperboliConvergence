
   template<typename T>
   void centralImplicitSmoothing( const Controls::GridControls1D&  grid,
                                  const Types::Real               alpha,
                                  const Array::Array1D<T>&           r0,
                                        Array::Array1D<T>&           r1 )
  {
      size_t i,n=grid.n;
      assert( r0.size() == n );
      assert( r1.size() == n );

      Array::Array1D<T>  wrk(n);

      Types::Real beta;
      int j;

      beta = 0.25*( alpha*alpha - 1. );

      bool fromZero=false;
      int m=100;

      r1 =r0;
      if( !fromZero ){ wrk=r0; }
      for( j=0; j<m; j++ )
     {
         r1=r0;

         for( i=0; i<n-1; i++ )
        {
            r1[i  ]+= beta*wrk[i+1];
            r1[i+1]+= beta*wrk[i  ];
        }
         for( i=1; i<n-1; i++ ){ r1[i]/=(1.+2.*beta); }

         if( grid.boundaryCondition=='p' )
        {
            r1[  0]+= beta*wrk[n-1];
            r1[n-1]+= beta*wrk[  0];

            r1[  0]/=(1.+2.*beta);
            r1[n-1]/=(1.+2.*beta);
        }
         else if( grid.boundaryCondition=='a' )
        {
            r1[  0]/=(1.+1.*beta);
            r1[n-1]/=(1.+1.*beta);
        }
         if( j<m-1 ){ wrk=r1; }
     }
      if( fromZero ){ r1=0.5*(r1+wrk); }

      r1*=alpha;
      return;
  }

