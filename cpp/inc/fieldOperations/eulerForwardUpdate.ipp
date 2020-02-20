# include <array1D/array1D.h>
# include <idealGas2D/idealGas2D.h>
# include <fieldOperations/fieldOperations.h>

   template<typename SolutionType>
   void eulerForwardUpdate( const IdealGas2D::Species& gas, float dt,
                            const Array1D< IdealGas2D::VariableSet< SolutionType > >& q0,
                                  Array1D< IdealGas2D::VariableSet< SolutionType > >& q,
                            const Array1D< IdealGas2D::ConservedDelta >&              r  )
  {
      int n=q0.size();

      assert( n= q.size() );
      assert( n= r.size() );

      IdealGas2D::ConservedVariables qc;

      IdealGas2D::VariableDelta< SolutionType > dq;

      for( int i=0; i<n; i++ )
     {
         qc = IdealGas2D::ConservedVariables( gas, q0[i] );
         qc+= dt*r[i];
         q[i] = IdealGas2D::VariableSet< SolutionType >( gas, qc );
     }
      return;
  }

