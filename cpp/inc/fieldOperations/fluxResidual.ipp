
   template<typename SolutionType, typename Flux, typename Limiter>
   void fluxResidual( const IdealGas2D::Species& gas, const char boundaryCondition,
                      const Flux flux, const Limiter limiter,
                      const Array::Array1D< IdealGas2D::VariableSet<   SolutionType > >&  q,
                      const Array::Array1D< IdealGas2D::VariableDelta< SolutionType > >& dq,
                            Array::Array1D< IdealGas2D::ConservedDelta >&                 r,
                                                             Types::Real&              lmax  )
  {
      int n=q.size();
      int i;

      assert( n = r.size() );

      Types::Real faceMetric[3]={1,0,1};
      Types::Real l,lm;

      IdealGas2D::ConservedDelta   f;
      IdealGas2D::VariableSet<SolutionType>     ql,     qr;
      IdealGas2D::VariableDelta<SolutionType>  dql,    dqr;
      IdealGas2D::VariableDelta<SolutionType>  dul,duc,dur;

      r = 0;
      l =-1;
      lm=-1;

      for( i=1; i<n; i++ )
     {
      // limit gradients
         dul = dq[i  ];
         duc = dq[i+1];
         dur = dq[i+2];

         limiter( dul,duc,dur, dql,dqr );

      // extrapolate face values
         ql = q[i-1] + 0.5*dql;
         qr = q[i  ] - 0.5*dqr;

      // evaluate flux
         flux( gas, faceMetric,
               IdealGas2D::State( gas, ql ),
               IdealGas2D::State( gas, qr ),
               f,l );

         lm  = fmax( lm, l );

      // accumulate residual
         r[i-1]-= f;
         r[i]  += f;
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
