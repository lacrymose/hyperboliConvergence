
namespace IdealGas2D
{
// default constructor
   template<typename VType>
   inline VariableDelta<VType>::VariableDelta()
  {
      var[0]=0.;
      var[1]=0.;
      var[2]=0.;
      var[3]=0.;
      return;
  }

// copy constructor
   template<typename VType>
   inline VariableDelta<VType>::VariableDelta(                     const VariableDelta<VType>& dq0 )
  {
      var[0]=dq0[0];
      var[1]=dq0[1];
      var[2]=dq0[2];
      var[3]=dq0[3];
      return;
  }

   // copy constructor to shadow linear transformation from another delta
   template<typename VType>
   inline VariableDelta<VType>::VariableDelta( const Species& gas, const State& state, const VariableDelta<VType>& dq0 )
  {
      var[0]=dq0[0];
      var[1]=dq0[1];
      var[2]=dq0[2];
      var[3]=dq0[3];
      return;
  }

// convert q -> dq
   template<typename VType>
   inline VariableDelta<VType>::VariableDelta(                     const VariableSet<VType>&    q0 )
  {
      var[0]=q0[0];
      var[1]=q0[1];
      var[2]=q0[2];
      var[3]=q0[3];
      return;
  }

// setter
   template<typename VType>
   inline VariableDelta<VType>& VariableDelta<VType>::operator=( Types::Real a )
  {
      var[0]=a;
      var[1]=a;
      var[2]=a;
      var[3]=a;
      return *this;
  }

// dq+=dq
   template<typename VType>
   inline VariableDelta<VType>& VariableDelta<VType>::operator+=( const VariableDelta<VType>& dq )
  {
      var[0]+=dq[0];
      var[1]+=dq[1];
      var[2]+=dq[2];
      var[3]+=dq[3];
      return *this;
  }

// dq-=dq
   template<typename VType>
   inline VariableDelta<VType>& VariableDelta<VType>::operator-=( const VariableDelta<VType>& dq )
  {
      var[0]-=dq[0];
      var[1]-=dq[1];
      var[2]-=dq[2];
      var[3]-=dq[3];
      return *this;
  }

// dq*=a
   template<typename VType>
   inline VariableDelta<VType>& VariableDelta<VType>::operator*=( Types::Real a )
  {
      var[0]*=a;
      var[1]*=a;
      var[2]*=a;
      var[3]*=a;
      return *this;
  }

// dq/=a
   template<typename VType>
   inline VariableDelta<VType>& VariableDelta<VType>::operator/=( Types::Real a )
  {
      Types::Real a1=1./a;
      var[0]*=a1;
      var[1]*=a1;
      var[2]*=a1;
      var[3]*=a1;
      return *this;
  }

// print dvariables to stream
   template<typename VType>
   inline std::ostream& operator<<( std::ostream& os, VariableDelta<VType> q )
  {
      os << std::left << std::setw(10) << q[0] << " "
         << std::left << std::setw(10) << q[1] << " "
         << std::left << std::setw(10) << q[2] << " "
         << std::left << std::setw(10) << q[3];
      return os;
  }
}

