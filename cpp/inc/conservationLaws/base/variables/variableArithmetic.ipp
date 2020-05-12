
// Arithmetic operations for VariableSets and VariableDeltas

   template<LawType Law, int nDim, BasisType<Law> Basis>    // dq=dq+dq
   inline VariableDelta<Law,nDim,Basis> operator+( const VariableDelta<Law,nDim,Basis>& dq0,
                                                   const VariableDelta<Law,nDim,Basis>& dq1 )
  {
      VariableDelta<Law,nDim,Basis> dq(dq0);
      dq+=dq1;
      return dq;
  }

   template<LawType Law, int nDim, BasisType<Law> Basis>    // dq=dq-dq
   inline VariableDelta<Law,nDim,Basis> operator-( const VariableDelta<Law,nDim,Basis>& dq0,
                                                   const VariableDelta<Law,nDim,Basis>& dq1 )
  {
      VariableDelta<Law,nDim,Basis> dq(dq0);
      dq-=dq1;
      return dq;
  }

   template<LawType Law, int nDim, BasisType<Law> Basis>     // dq=q-q
   inline VariableDelta<Law,nDim,Basis> operator-( const VariableSet<Law,nDim,Basis>& q0,
                                                   const VariableSet<Law,nDim,Basis>& q1 )
  {
      using Delta = VariableDelta<Law,nDim,Basis>;
      Delta dq = Delta(q0) - Delta(q1);
      return dq;
  }

   template<LawType Law, int nDim, BasisType<Law> Basis>     //  q=q+dq
   inline VariableSet<Law,nDim,Basis>   operator+( const VariableSet<  Law,nDim,Basis>&  q0,
                                                   const VariableDelta<Law,nDim,Basis>& dq0 )
  {
      VariableSet<Law,nDim,Basis> q(q0);
      q+=dq0;
      return q;
  }

   template<LawType Law, int nDim, BasisType<Law> Basis>     //  q=q-dq
   inline VariableSet<Law,nDim,Basis>   operator-( const VariableSet<  Law,nDim,Basis>&  q0,
                                                   const VariableDelta<Law,nDim,Basis>& dq0 )
  {
      VariableSet<Law,nDim,Basis> q(q0);
      q-=dq0;
      return q;
  }

   template<LawType Law, int nDim, BasisType<Law> Basis>    //  dq=a*dq
   inline VariableDelta<Law,nDim,Basis> operator*( const Types::Real a,
                                                   const VariableDelta<Law,nDim,Basis>& dq0 )
  {
      VariableDelta<Law,nDim,Basis> dq(dq0);
      dq*=a;
      return dq;
  }

   template<LawType Law, int nDim, BasisType<Law> Basis>    // dq=dq*a
   inline VariableDelta<Law,nDim,Basis> operator*( const VariableDelta<Law,nDim,Basis>& dq0,
                                                   const Types::Real a )
  {
      VariableDelta<Law,nDim,Basis> dq(dq0);
      dq*=a;
      return dq;
  }

   template<LawType Law, int nDim, BasisType<Law> Basis>    // dq=dq/a
   inline VariableDelta<Law,nDim,Basis> operator/( const VariableDelta<Law,nDim,Basis>& dq0,
                                                   const Types::Real a )
  {
      VariableDelta<Law,nDim,Basis> dq(dq0);
      dq/=a;
      return dq;
  }

