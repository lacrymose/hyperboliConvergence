
namespace IdealGas2D
{
// Arithmetic operations for VariableSets and VariableDeltas

   template<typename VType>    // dq=dq+dq
   inline VariableDelta<VType> operator+( const VariableDelta<VType>& dq0,
                                          const VariableDelta<VType>& dq1 )
  {
      VariableDelta<VType> dq2(dq0);
      dq2+=dq1;
      return dq2;
  }

   template<typename VType>    // dq=dq-dq
   inline VariableDelta<VType> operator-( const VariableDelta<VType>& dq0,
                                          const VariableDelta<VType>& dq1 )
  {
      VariableDelta<VType> dq2(dq0);
      dq2-=dq1;
      return dq2;
  }

   template<typename VType>     // dq=q-q
   inline VariableDelta<VType> operator-( const VariableSet<VType>& q0,
                                          const VariableSet<VType>& q1 )
  {
      VariableDelta<VType> dq(q0);
      dq-=q1;
      return dq;
  }

   template<typename VType>     //  q=q+dq
   inline VariableSet<VType>   operator+( const VariableSet<VType>&    q0,
                                          const VariableDelta<VType>& dq0 )
  {
      VariableSet<VType> q1(q0);
      q1+=dq0;
      return q1;
  }

   template<typename VType>     //  q=q-dq
   inline VariableSet<VType>   operator-( const VariableSet<VType>&    q0,
                                          const VariableDelta<VType>& dq0 )
  {
      VariableSet<VType> q1(q0);
      q1-=dq0;
      return q1;
  }

   template<typename VType>    //  dq=a*dq
   inline VariableDelta<VType> operator*(                      float    a,
                                          const VariableDelta<VType>& dq0 )
  {
      VariableDelta<VType> dq1(dq0);
      dq1*=a;
      return dq1;
  }

   template<typename VType>    // dq=dq*a
   inline VariableDelta<VType> operator*( const VariableDelta<VType>& dq0,
                                                               float    a )
  {
      VariableDelta<VType> dq1(dq0);
      dq1*=a;
      return dq1;
  }

   template<typename VType>    // dq=dq/a
   inline VariableDelta<VType> operator/( const VariableDelta<VType>& dq0,
                                                               float    a )
  {
      VariableDelta<VType> dq1(dq0);
      dq1/=a;
      return dq1;
  }

// Arithmetic operations for arrays of VariableSets and VariableDeltas

//

// dq=q-q
   template<typename VType>
   inline Array::Array1D<VariableDelta<VType>> operator-( const Array::Array1D<VariableSet<VType>>& q0,
                                                          const Array::Array1D<VariableSet<VType>>& q1 )
  {
      int n=q0.size();
      assert( n==q1.size() );
      Array::Array1D<VariableDelta<VType>> dq(n);
      for( int i=0; i<n; i++ )
     {
         dq[i]=q0[i]-q1[i];
     }
      return dq;
  }

// q=q+dq
   template<typename VType>
   inline Array::Array1D<VariableSet<VType>> operator+( const Array::Array1D<VariableSet<  VType>>& q0,
                                                        const Array::Array1D<VariableDelta<VType>>& dq )
  {
      int n=q0.size();
      assert( n==dq.size() );
      Array::Array1D<VariableSet<VType>> q1(n);
      q1 =q0;
      q1+=dq;
      return q1;
  }

// q=q-dq
   template<typename VType>
   inline Array::Array1D<VariableSet<VType>> operator-( const Array::Array1D<VariableSet<  VType>>& q0,
                                                        const Array::Array1D<VariableDelta<VType>>& dq )
  {
      int n=q0.size();
      assert( n==dq.size() );
      Array::Array1D<VariableSet<VType>> q1(n);
      q1 =q0;
      q1-=dq;
      return q1;
  }

// q+=dq
   template<typename VType>
   inline Array::Array1D<VariableSet<VType>>& operator+=(       Array::Array1D<VariableSet<  VType>>& q,
                                                          const Array::Array1D<VariableDelta<VType>>& dq )
  {
      int n=q.size();
      assert( n==dq.size() );
      for( int i=0; i<n; i++ )
     {
         q[i]+=dq[i];
     }
      return q;
  }

// q-=dq
   template<typename VType>
   inline Array::Array1D<VariableSet<VType>>& operator-=(       Array::Array1D<VariableSet<  VType>>& q,
                                                          const Array::Array1D<VariableDelta<VType>>& dq )
  {
      int n=q.size();
      assert( n==dq.size() );
      for( int i=0; i<n; i++ )
     {
         q[i]-=dq[i];
     }
      return q;
  }
}
