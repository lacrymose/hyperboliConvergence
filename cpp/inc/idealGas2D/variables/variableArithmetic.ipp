
namespace IdealGas2D
{
// Arithmetic operations for VariableSets and VariableDeltas
   template<typename VType>     // dq=q-q
   inline VariableDelta<VType> operator-( const VariableSet<VType>& q0, const VariableSet<VType>& q1 )
  {
      VariableDelta<VType> dq(q0);
      dq-=q1;
      return dq;
  }

   template<typename VType>    // dq=dq+dq
   inline VariableDelta<VType> operator+( const VariableDelta<VType>& dq0, const VariableDelta<VType>& dq1 )
  {
      VariableDelta<VType> dq2(dq0);
      dq2+=dq1;
      return dq2;
  }

   template<typename VType>    // dq=dq-dq
   inline VariableDelta<VType> operator-( const VariableDelta<VType>& dq0, const VariableDelta<VType>& dq1 )
  {
      VariableDelta<VType> dq2(dq0);
      dq2-=dq1;
      return dq2;
  }

   template<typename VType>     //  q=q+dq
   inline VariableSet<VType>   operator+( const VariableSet<VType>&    q0, const VariableDelta<VType>& dq0 )
  {
      VariableSet<VType> q1(q0);
      q1+=dq0;
      return q1;
  }

   template<typename VType>     //  q=q-dq
   inline VariableSet<VType>   operator-( const VariableSet<VType>&    q0, const VariableDelta<VType>& dq0 )
  {
      VariableSet<VType> q1(q0);
      q1-=dq0;
      return q1;
  }

   template<typename VType>    //  dq=a*dq
   inline VariableDelta<VType> operator*(                      float    a, const VariableDelta<VType>& dq0 )
  {
      VariableDelta<VType> dq1(dq0);
      dq1*=a;
      return dq1;
  }

   template<typename VType>    // dq=dq*a
   inline VariableDelta<VType> operator*( const VariableDelta<VType>& dq0,                      float    a )
  {
      VariableDelta<VType> dq1(dq0);
      dq1*=a;
      return dq1;
  }

   template<typename VType>    // dq=dq/a
   inline VariableDelta<VType> operator/( const VariableDelta<VType>& dq0,                      float    a )
  {
      VariableDelta<VType> dq1(dq0);
      dq1/=a;
      return dq1;
  }

}
