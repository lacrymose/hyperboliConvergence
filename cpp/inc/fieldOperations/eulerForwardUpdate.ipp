# include <array1D/array1D.h>
# include <idealGas2D/idealGas2D.h>
# include <fieldOperations/fieldOperations.h>

   template<typename SolutionType>
   void eulerForwardUpdate( const IdealGas2D::Species& gas, float dt,
                            const Array::Array1D< IdealGas2D::VariableSet< SolutionType > >& q0,
                                  Array::Array1D< IdealGas2D::VariableSet< SolutionType > >& q,
                            const Array::Array1D< IdealGas2D::ConservedDelta >&              r,
                                                    IdealGas2D::ConservedDelta&              res  )
  {
      int n=q0.size();
      int i,j;

      assert( n= q.size() );
      assert( n= r.size() );

      IdealGas2D::ConservedDelta               f;
      IdealGas2D::ConservedVariables          qc;
      IdealGas2D::VariableDelta<SolutionType> dq;

      res=0;
      for( i=0; i<n; i++ )
     {
         qc   = IdealGas2D::ConservedVariables( gas, q0[i] );
         qc  += dt*r[i];
         q[i] = IdealGas2D::VariableSet< SolutionType >( gas, qc );

//       dq = IdealGas2D::VariableDelta<SolutionType>( gas, IdealGas2D::State(gas,q[i]), r[i] );
//       q[i]+=dt*dq;

         for( j=0; j<4; j++ ){ res[j] += fabs( dt*r[i][j] ); }
     }
      return;
  }

