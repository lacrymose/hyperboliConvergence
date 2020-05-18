

# include <types.h>

# include <ostream>

// p+=dir
   template<int nDim>
   Direction<nDim>& Direction<nDim>::operator+=( const Direction& d )
  {
      for( int i=0; i<nDim; i++ ){ x[i]+=d.x[i]; }
      return *this;
  }

// p-=dir
   template<int nDim>
   Direction<nDim>& Direction<nDim>::operator-=( const Direction& d )
  {
      for( int i=0; i<nDim; i++ ){ x[i]-=d.x[i]; }
      return *this;
  }

// p=a
   template<int nDim>
   Direction<nDim>& Direction<nDim>::operator =( const Types::Real a )
  {
      for( int i=0; i<nDim; i++ ){ x[i]=a; }
      return *this;
  }

// p*=a
   template<int nDim>
   Direction<nDim>& Direction<nDim>::operator*=( const Types::Real a )
  {
      for( int i=0; i<nDim; i++ ){ x[i]*=a; }
      return *this;
  }

// p/=a
   template<int nDim>
   Direction<nDim>& Direction<nDim>::operator/=( const Types::Real a )
  {
      const Types::Real a1=1./a;
      for( int i=0; i<nDim; i++ ){ x[i]*=a1; }
      return *this;
  }

   template<int nDim>
   std::ostream& operator<<( std::ostream& os, const Direction<nDim>& d )
  {
      for( int i=0; i<nDim-1; i++ )
     {
         os << d[i] << " ";
     }
      os << d[nDim-1];
      return os;
  }

