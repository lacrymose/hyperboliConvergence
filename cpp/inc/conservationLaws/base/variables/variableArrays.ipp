
# include <array1D/array1D.h>

// Arithmetic operations for arrays of VariableSets and VariableDeltas

// dq=q-q
   template<LawType Law, int nDim, BasisType<Law> Basis>
   inline Array::Array1D<VariableDelta<Law,nDim,Basis>> operator-( const Array::Array1D<VariableSet<Law,nDim,Basis>>& q0,
                                                                   const Array::Array1D<VariableSet<Law,nDim,Basis>>& q1 )
  {
      size_t n=q0.size();
      assert( n==q1.size() );

      Array::Array1D<VariableDelta<Law,nDim,Basis>> dq(n);
      for( size_t i=0; i<n; i++ )
     {
         dq[i]=q0[i]-q1[i];
     }
      return dq;
  }

// q=q+dq
   template<LawType Law, int nDim, BasisType<Law> Basis>
   inline Array::Array1D<VariableSet<Law,nDim,Basis>> operator+( const Array::Array1D<VariableSet<  Law,nDim,Basis>>& q0,
                                                                 const Array::Array1D<VariableDelta<Law,nDim,Basis>>& dq )
  {
      size_t n=q0.size();
      assert( n==dq.size() );

      Array::Array1D<VariableSet<Law,nDim,Basis>> q1(q0);
      for( size_t i=0; i<n; i++ )
     {
         q1[i]+=dq[i];
     }
      return q1;
  }

// q=q-dq
   template<LawType Law, int nDim, BasisType<Law> Basis>
   inline Array::Array1D<VariableSet<Law,nDim,Basis>> operator-( const Array::Array1D<VariableSet<  Law,nDim,Basis>>& q0,
                                                                 const Array::Array1D<VariableDelta<Law,nDim,Basis>>& dq )
  {
      size_t n=q0.size();
      assert( n==dq.size() );

      Array::Array1D<VariableSet<Law,nDim,Basis>> q1(q0);
      for( size_t i=0; i<n; i++ )
     {
         q1[i]-=dq[i];
     }
      return q1;
  }

// q+=dq
   template<LawType Law, int nDim, BasisType<Law> Basis>
   inline Array::Array1D<VariableSet<Law,nDim,Basis>>& operator+=(       Array::Array1D<VariableSet<  Law,nDim,Basis>>& q,
                                                                   const Array::Array1D<VariableDelta<Law,nDim,Basis>>& dq )
  {
      size_t n=q.size();
      assert( n==dq.size() );

      for( size_t i=0; i<n; i++ )
     {
         q[i]+=dq[i];
     }
      return q;
  }

// q-=dq
   template<LawType Law, int nDim, BasisType<Law> Basis>
   inline Array::Array1D<VariableSet<Law,nDim,Basis>>& operator-=(       Array::Array1D<VariableSet<  Law,nDim,Basis>>& q,
                                                                   const Array::Array1D<VariableDelta<Law,nDim,Basis>>& dq )
  {
      size_t n=q.size();
      assert( n==dq.size() );

      for( size_t i=0; i<n; i++ )
     {
         q[i]-=dq[i];
     }
      return q;
  }

