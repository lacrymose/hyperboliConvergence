
namespace IdealGas2D
{
// default constructor
   template<typename VType>
   inline VariableSet<VType>::VariableSet()
  {
      var[0]=0.;
      var[1]=0.;
      var[2]=0.;
      var[3]=0.;
      return;
  }

// copy constructor
   template<typename VType>
   inline VariableSet<VType>::VariableSet(                     const VariableSet<VType>& q0 )
  {
      var[0]=q0[0];
      var[1]=q0[1];
      var[2]=q0[2];
      var[3]=q0[3];
      return;
  }

// copy constructor
   template<typename VType>
   inline VariableSet<VType>::VariableSet( const Species& gas, const VariableSet<VType>& q0 )
  {
      var[0]=q0[0];
      var[1]=q0[1];
      var[2]=q0[2];
      var[3]=q0[3];
      return;
  }

// convert dq -> q
   template<typename VType>
   inline VariableSet<VType>::VariableSet(                     const VariableDelta<VType>& dq0 )
  {
      var[0]=dq0[0];
      var[1]=dq0[1];
      var[2]=dq0[2];
      var[3]=dq0[3];
      return;
  }

// nonlinear transformation from VType2 to VType
   template<typename VType>
   template<typename VType2>
   inline VariableSet<VType>::VariableSet( const Species& gas, const VariableSet<VType2>& q0 )
  {
      static_assert( CheckTypes<VType,Conserved>::val ||
                     CheckTypes<VType,Viscous  >::val,
                    "\n\nWarning:\n"
                    "VariableSet<VType>::VariableSet( const Species& gas, const VariableSet<VType2>& q0 )\n"
                    "is not yet defined for these VariableTypes\n" );
      assert( false );
  }

// nonlinear transformation from State to VType
   template<typename VType>
   inline VariableSet<VType>::VariableSet( const Species& gas, const State& state )
  {
      static_assert( CheckTypes<VType,Conserved>::val ||
                     CheckTypes<VType,Viscous  >::val,
                    "\n\nWarning:\n"
                    "VariableSet<VType>::VariableSet( const Species& gas, const State& s0 )\n"
                    "is not yet defined for these VariableTypes\n" );
      assert( false );
  }

// setter
   template<typename VType>
   inline VariableSet<VType>& VariableSet<VType>::operator=( float a )
  {
      var[0]=a;
      var[1]=a;
      var[2]=a;
      var[3]=a;
      return *this;
  }

// q+=dq
   template<typename VType>
   inline VariableSet<VType>& VariableSet<VType>::operator+=( const VariableDelta<VType>& dq )
  {
      var[0]+=dq[0];
      var[1]+=dq[1];
      var[2]+=dq[2];
      var[3]+=dq[3];
      return *this;
  }

// q-=dq
   template<typename VType>
   inline VariableSet<VType>& VariableSet<VType>::operator-=( const VariableDelta<VType>& dq )
  {
      var[0]-=dq[0];
      var[1]-=dq[1];
      var[2]-=dq[2];
      var[3]-=dq[3];
      return *this;
  }

// print variables to stream
   template<typename VType>
   inline std::ostream& operator<<( std::ostream& os, VariableSet<VType> q )
  {
      os << q[0] << " "
         << q[1] << " "
         << q[2] << " "
         << q[3];
      return os;
  }
}

