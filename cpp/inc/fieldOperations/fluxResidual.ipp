
//typedef IdealGas2D::Ausm            Flux;
//typedef IdealGas2D::LaxFriedrichs   Flux;

   template<char C, typename Flux>
   void fluxResidual( const IdealGas2D::Species& gas, const char boundaryCondition, Flux flux,
                            const Array1D< IdealGas2D::VariableSet<C> >&       q,
                                  Array1D< IdealGas2D::ConservedVariables >&   r,
                                  float& lmax )
  {
      int n=q.size();

      assert( n = r.size() );

      float faceMetric[3]={1,0,1};
      float l,lm;

      IdealGas2D::ConservedVariables   f;

      r =0.;
      l =-1;
      lm=-1;

      for( int i=0; i<n-1; i++ )
     {
         flux( gas, faceMetric,
               IdealGas2D::State( gas, q[i]   ),
               IdealGas2D::State( gas, q[i+1] ),
               f,l );

         lm = fmax( lm, l );

         r[i]  -= f;
         r[i+1]+= f;
     }

   // zero flux (adiabatic) boundaries
      if( boundaryCondition=='a' )
     {
         r[0]=0.;
         r[n-1]=0.;
     }
   // periodic boundaries
      else if( boundaryCondition=='p' )
     {
         flux( gas, faceMetric,
               IdealGas2D::State( gas, q[n-1] ),
               IdealGas2D::State( gas, q[0]   ),
               f,l );

         lm = fmax( lm, l );

         r[n-1]-= f;
         r[0]  += f;
     }

      lmax = lm;

      return;
  }
