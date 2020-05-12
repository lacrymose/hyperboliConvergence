# ifndef GEOMETRY_H
# define GEOMETRY_H

# include <types.h>

# include <array>

namespace Geometry
{
/*
 * {Point,Direction} form an affine space
 */
   template<int nDim>
   struct Point
  {
      std::array<Types::Real,nDim>  x;

            Types::Real& operator[](       int i );
      const Types::Real& operator[]( const int i ) const;

      Point& operator+=( const Direction&  d );
      Point& operator-=( const Direction&  d );
      Point& operator =( const Types::Real a );
  };

   template<int nDim>
   struct Direction
  {
      std::array<Types::Real,nDim>  x;

            Types::Real& operator[](       int i );
      const Types::Real& operator[]( const int i ) const;

      Direction& operator+=( const Direction&  d );
      Direction& operator-=( const Direction&  d );
      Direction& operator*=( const Types::Real a );
      Direction& operator/=( const Types::Real a );
      Direction& operator =( const Types::Real a );
  };

/*
 * Metric for a surface in nDim dimensions. Direction m[0] is surface normal, m[1...] are surface tangents
 */
   template<int nDim>
   struct Metric
  {
      std::array<Direction<nDim>,nDim>  m;

            Direction<nDim>& operator[](       int i );
      const Direction<nDim>& operator[]( const int i ) const;
  };

// used to represent sections of a space of the same # of dimensions as the space
   template<int nDim>
   struct Volume
  {
      Types::Real  volume;
      Point<nDim>  centre;
  };

// used to represent hypersurfaces in a space
   template<int nDim>
   struct Surface
  {
      Types::Real     area;
      Point<nDim>   centre;
      Metric<nDim>  metric;
  };

   Surface<2> surface( std::array<Point,2> );
   Surface<3> surface( std::array<Point,4> );

   Volume<2> volume( std::array<Point,4> );
   Volume<3> volume( std::array<Point,8> );
}

# endif
