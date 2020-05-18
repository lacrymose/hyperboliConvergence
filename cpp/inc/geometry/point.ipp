
# include <types.h>

# include <ostream>

// p+=dir
   template<int nDim>
   Point<nDim>& Point<nDim>::operator+=( const Direction& d )
  {
      for( int i=0; i<nDim; i++ ){ x[i]+=d.x[i]; }
      return *this;
  }

// p-=dir
   template<int nDim>
   Point<nDim>& Point<nDim>::operator-=( const Direction& d )
  {
      for( int i=0; i<nDim; i++ ){ x[i]-=d.x[i]; }
      return *this;
  }

// p=a
   template<int nDim>
   Point<nDim>& Point<nDim>::operator =( const Types::Real a )
  {
      for( int i=0; i<nDim; i++ ){ x[i]=a; }
      return *this;
  }

   template<int nDim>
   std::ostream& operator<<( std::ostream& os, const Point<nDim>& p )
  {
      for( int i=0; i<nDim-1; i++ )
     {
         os << p[i] << " ";
     }
      os << p[nDim-1];
      return os;
  }

