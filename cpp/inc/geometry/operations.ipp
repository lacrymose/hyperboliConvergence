
   template<int nDim>
   Types::Real length2( const Direction<nDim>& d )
  {
      Types::Real l=0;
      for( const Types::Real& x : d.x ){ l+=x*x; }
      return l;
  }

   template<int nDim>
   Types::Real length( const Direction<nDim>& d )
  {
      return std::sqrt(length2(d));
  }

   Direction<2> cross( const Direction<2>& a )
  {
      return Direction<2>{-a[1],a[0]};
  }

   Direction<3> cross( const Direction<3>& a,
                       const Direction<3>& b )
  {
      return Direction<3>{ a[1]*b[2] - a[2]*b[1],
                           a[2]*b[0] - a[0]*b[2],
                           a[0]*b[1] - a[1]*b[0] };
  }

