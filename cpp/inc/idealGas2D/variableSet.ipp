
# include <iostream>
# include <assert.h>

namespace IdealGas2D
{
   template< char C >
   inline VariableSet<C>::VariableSet()
  {
      var[0]=0.;
      var[1]=0.;
      var[2]=0.;
      var[3]=0.;
      return;
  }

   template< char C >
   inline VariableSet<C>::VariableSet( const Species& gas, const VariableSet<C>& q0 )
  {
      var[0]=q0[0];
      var[1]=q0[1];
      var[2]=q0[2];
      var[3]=q0[3];
      return;
  }

   template< char C >
   template< char D >
   inline VariableSet<C>::VariableSet( const Species& gas, const VariableSet<D>& q0 )
  {
      std::cout << std::endl;
      std::cout << "Warning:" << std::endl;
      std::cout << "VariableSet<C>::VariableSet( const VariableSet<D>& q0 )"
                << " is not yet defined for "
                << "C = '" << C << "' and "
                << "D = '" << D << "'" << std::endl;
      std::cout << std::endl;
      assert( false );
  }

   template< char C >
   inline VariableSet<C>::VariableSet( const Species& gas, const State& s0 )
  {
      std::cout << std::endl;
      std::cout << "Warning:" << std::endl;
      std::cout << "VariableSet<C>::VariableSet( const State& s0 )"
                << " is not yet defined for "
                << "C = '" << C << "'" << std::endl;
      std::cout << std::endl;
      assert( false );
  }

   template< char C >
   inline VariableSet<C>& VariableSet<C>::operator=( float d )
  {
      var[0]=d;
      var[1]=d;
      var[2]=d;
      var[3]=d;
      return *this;
  }

   template< char C >
   inline VariableSet<C>& VariableSet<C>::operator+=( const VariableSet<C>& q )
  {
      var[0]+=q[0];
      var[1]+=q[1];
      var[2]+=q[2];
      var[3]+=q[3];
      return *this;
  }

   template< char C >
   inline VariableSet<C>& VariableSet<C>::operator-=( const VariableSet<C>& q )
  {
      var[0]-=q[0];
      var[1]-=q[1];
      var[2]-=q[2];
      var[3]-=q[3];
      return *this;
  }

   template< char C >
   inline VariableSet<C>& VariableSet<C>::operator*=( float d )
  {
      var[0]*=d;
      var[1]*=d;
      var[2]*=d;
      var[3]*=d;
      return *this;
  }

   template< char C >
   inline VariableSet<C>& VariableSet<C>::operator/=( float d )
  {
      float d1=1./d;
      var[0]*=d1;
      var[1]*=d1;
      var[2]*=d1;
      var[3]*=d1;
      return *this;
  }

   template< char C >
   inline VariableSet<C> operator+( const VariableSet<C>& q0, const VariableSet<C>& q1 )
  {
      VariableSet<C> q2;
      q2[0] = q0[0] + q1[0];
      q2[1] = q0[1] + q1[1];
      q2[2] = q0[2] + q1[2];
      q2[3] = q0[3] + q1[3];
      return q2;
  }

   template< char C >
   inline VariableSet<C> operator-( const VariableSet<C>& q0, const VariableSet<C>& q1 )
  {
      VariableSet<C> q2;
      q2[0] = q0[0] - q1[0];
      q2[1] = q0[1] - q1[1];
      q2[2] = q0[2] - q1[2];
      q2[3] = q0[3] - q1[3];
      return q2;
  }

   template< char C >
   inline VariableSet<C> operator*( const VariableSet<C>& q0, float d )
  {
      VariableSet<C> q1;
      q1[0] = d*q0[0];
      q1[1] = d*q0[1];
      q1[2] = d*q0[2];
      q1[3] = d*q0[3];
      return q1;
  }

   template< char C >
   inline VariableSet<C> operator*( float d, const VariableSet<C>& q0 )
  {
      VariableSet<C> q1;
      q1[0] = d*q0[0];
      q1[1] = d*q0[1];
      q1[2] = d*q0[2];
      q1[3] = d*q0[3];
      return q1;
  }

   template< char C >
   inline VariableSet<C> operator/( const VariableSet<C>& q0, float d )
  {
      VariableSet<C> q1;
      float d1=1./d;
      q1[0] = q0[0]*d1;
      q1[1] = q0[1]*d1;
      q1[2] = q0[2]*d1;
      q1[3] = q0[3]*d1;
      return q1;
  }

   template< char C >
   inline std::ostream& operator<<( std::ostream& os, VariableSet<C> q )
  {
      os << q[0] << " "
         << q[1] << " "
         << q[2] << " "
         << q[3];
      return os;
  }
}

