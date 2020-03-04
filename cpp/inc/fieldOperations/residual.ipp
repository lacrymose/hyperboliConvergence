
   template<typename SolutionType, typename ResidualType>
   void residual( const IdealGas2D::Species& gas,
                  const Array::Array1D< IdealGas2D::VariableSet<  SolutionType > >&   q0,
                  const Array::Array1D< IdealGas2D::VariableSet<  SolutionType > >&   q1,
                                        IdealGas2D::VariableDelta<ResidualType>&     res )
  {
      assert( q0.size() == q1.size() );

      IdealGas2D::VariableDelta<ResidualType> dq;

      res=0;
      for( size_t i=0; i<q0.size(); i++ )
     {
         dq =  IdealGas2D::VariableSet<ResidualType>( gas, q0[i] )
             - IdealGas2D::VariableSet<ResidualType>( gas, q1[i] );

         res+= IdealGas2D::abs( dq );
     }
  }
