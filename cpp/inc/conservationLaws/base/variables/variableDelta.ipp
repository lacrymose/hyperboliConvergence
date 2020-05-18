
# include <types.h>

# include <ostream>


//// explicit copy construction from VariableDelta
//   template<LawType Law, int nDim, BasisType<Law> Basis>
//   explicit VariableDelta<Law,nDim,Basis>::VariableDelta( const VariableSet& q ) noexcept
//  {
//      for( int i=0; i<(nVar<Law,nDim>); i++ ){ var[i]=q[i]; }
//  }
//
//// explicit copy construction from VariableDelta
//   template<LawType Law, int nDim, BasisType<Law> Basis>
//   explicit VariableDelta<Law,nDim,Basis>::VariableDelta( VariableSet&& q ) noexcept
//  {
//      var( std::move(q.var) );
//  }

// in-place addition
   template<LawType Law, int nDim, BasisType<Law> Basis>
   VariableDelta<Law,nDim,Basis>& VariableDelta<Law,nDim,Basis>::operator+=( const VariableDelta& dq0 )
  {
      for( int i=0; i<(nVar<Law,nDim>); i++ ){ var[i]+=dq0[i]; }
      return *this;
  }

// in-place subtraction
   template<LawType Law, int nDim, BasisType<Law> Basis>
   VariableDelta<Law,nDim,Basis>& VariableDelta<Law,nDim,Basis>::operator-=( const VariableDelta& dq0 )
  {
      for( int i=0; i<(nVar<Law,nDim>); i++ ){ var[i]-=dq0[i]; }
      return *this;
  }

// in-place multiplication
   template<LawType Law, int nDim, BasisType<Law> Basis>
   VariableDelta<Law,nDim,Basis>& VariableDelta<Law,nDim,Basis>::operator*=( const Types::Real a )
  {
      for( int i=0; i<(nVar<Law,nDim>); i++ ){ var[i]*=a; }
      return *this;
  }

// in-place division
   template<LawType Law, int nDim, BasisType<Law> Basis>
   VariableDelta<Law,nDim,Basis>& VariableDelta<Law,nDim,Basis>::operator/=( const Types::Real a )
  {
      const Types::Real a1=1./a;
      for( int i=0; i<(nVar<Law,nDim>); i++ ){ var[i]*=a1; }
      return *this;
  }

// assignment
   template<LawType Law, int nDim, BasisType<Law> Basis>
   VariableDelta<Law,nDim,Basis>& VariableDelta<Law,nDim,Basis>::operator =( const Types::Real a )
  {
      for( int i=0; i<(nVar<Law,nDim>); i++ ){ var[i]=a; }
      return *this;
  }

// print variable deltas to stream
   template<LawType Law, int nDim, BasisType<Law> Basis>
   std::ostream& operator<<( std::ostream& os, const VariableDelta<Law,nDim,Basis> dq )
  {
      for( int i=0; i<(nVar<Law,nDim>)-1; i++ )
     {
         os << dq[i] << " ";
     }
      os << dq[(nVar<Law,nDim>)-1];
      return os;
  }

