
   template<typename T>
   void multigridRestriction( const Controls::GridControls1D&  grid,
                              const Array::Array1D<T>&           rh,
                                    Array::Array1D<T>&          r2h )
  {
      size_t i,n=grid.n;
      assert( rh.size( ) == n   );
      assert( r2h.size() == n/2 );

      for( i=0; i<n/2; i++ )
     {
      // agglomeration
         r2h[i] = rh[2*i  ];
         r2h[i]+= rh[2*i+1];
     }

      return;
  }

   template<typename T>
   void multigridProlongation( const Controls::GridControls1D&  grid,
                               const Array::Array1D<T>&          e2h,
                                     Array::Array1D<T>&           eh )
  {
      size_t i,n=grid.n;
      assert( eh.size( ) == n   );
      assert( e2h.size() == n/2 );

      auto de = eh[1]-eh[0];

      for( i=1; i<n/2-1; i++ )
     {
      // central difference gradient
         de = 0.5*(e2h[i+1] - e2h[i-1]);

      // linear interpolation
         eh[2*i  ] = e2h[i] - 0.25*de;
         eh[2*i+1] = e2h[i] + 0.25*de;
     }

      if( grid.boundaryCondition=='p' )
     {
         de = 0.5*(e2h[1] - e2h[n/2-1]);

         eh[0] = e2h[0] - 0.25*de;
         eh[1] = e2h[0] + 0.25*de;

         de = 0.5*(e2h[0] - e2h[n/2-2]);

         eh[n-2] = e2h[n/2-1] - 0.25*de;
         eh[n-1] = e2h[n/2-1] + 0.25*de;
     }
      else if( grid.boundaryCondition=='a' )
     {
         de = e2h[1]-e2h[0];

         eh[0] = e2h[0] - 0.25*de;
         eh[1] = e2h[0] + 0.25*de;

         de = e2h[n/2-1] - e2h[n/2-2];

         eh[n-2] = e2h[n/2-1] - 0.25*de;
         eh[n-1] = e2h[n/2-1] + 0.25*de;
     }

      return;
  }

