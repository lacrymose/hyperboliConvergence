
# include <types.h>

# include <ostream>

//// explicit copy construction from VariableDelta
//   template<LawType Law, int nDim, BasisType<Law> Basis>
//   explicit inline VariableSet<Law,nDim,Basis>::VariableSet( const VariableDelta& dq ) noexcept
//  {
//      for( int i=0; i<(nVar<Law,nDim>); i++ ){ var[i]=dq[i]; }
//  }
//
//// explicit copy construction from VariableDelta
//   template<LawType Law, int nDim, BasisType<Law> Basis>
//   explicit inline VariableSet<Law,nDim,Basis>::VariableSet( VariableDelta&& dq ) noexcept
//  {
//      var( std::move(dq.var) );
//  }

// q+=dq
   template<LawType Law, int nDim, BasisType<Law> Basis>
   VariableSet<Law,nDim,Basis>& VariableSet<Law,nDim,Basis>::operator+=( const VariableDelta& dq0 )
  {
      for( int i=0; i<(nVar<Law,nDim>); i++ ){ var[i]+=dq0[i]; }
      return *this;
  }

// q-=dq
   template<LawType Law, int nDim, BasisType<Law> Basis>
   VariableSet<Law,nDim,Basis>& VariableSet<Law,nDim,Basis>::operator-=( const VariableDelta& dq0 )
  {
      for( int i=0; i<(nVar<Law,nDim>); i++ ){ var[i]-=dq0[i]; }
      return *this;
  }

// q=a
   template<LawType Law, int nDim, BasisType<Law> Basis>
   VariableSet<Law,nDim,Basis>& VariableSet<Law,nDim,Basis>::operator =( const Types::Real a )
  {
      for( int i=0; i<(nVar<Law,nDim>); i++ ){ var[i]=a; }
      return *this;
  }

// print variable set to stream
   template<LawType Law, int nDim, BasisType<Law> Basis>
   std::ostream& operator<<( std::ostream& os, const VariableSet<Law,nDim,Basis>& q )
  {
      for( int i=0; i<(nVar<Law,nDim>)-1; i++ )
     {
         os << q[i] << " ";
     }
      os << q[(nVar<Law,nDim>)-1];
      return os;
  }

