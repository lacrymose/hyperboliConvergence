
# pragma once

# include <utils/concepts.h>
# include <types.h>

# include <array>
# include <vector>
# include <ostream>

# include <cmath>

namespace Geometry
{
// declarations
   template<int nDim> struct Point;
   template<int nDim> struct Direction;
   template<int nDim> struct Metric;
   template<int nDim> struct Volume;
   template<int nDim> struct Surface;

/*
 * {Point,Direction} form an affine space
 */
   template<int nDim>
   struct Point
  {
      using Direction = Geometry::Direction<nDim>;

      std::array<Types::Real,nDim>  x{};

   // default, copy and move constructors
      Point() = default;
      Point( const Point&  ) = default;
      Point(       Point&& ) = default;

   // only explicit construction from VariableDelta
      explicit Point( const Direction&  d ) noexcept : x(d.x) {}
      explicit Point(       Direction&& d ) noexcept : x(std::move(d.x)) {}

   // copy/move assignment
      Point& operator=( const Point&  ) = default;
      Point& operator=(       Point&& ) = default;

   // initialiser list constructor, must be a length nDim list of Types::Real
      template<typename... T>
         requires   Same<Types::Real,T...>
                 && is_pack_of_n<nDim,T...>::value
      Point( T... r ) noexcept : x{r...} {}

   // accessors
            Types::Real& operator[]( const int i )       { return x[i]; }
      const Types::Real& operator[]( const int i ) const { return x[i]; }

   // in-place arithmetic
      Point& operator+=( const Direction&  d );
      Point& operator-=( const Direction&  d );
      Point& operator =( const Types::Real a );
  };

/*
 * send each element of Point to stream, seperated by a single space
 */
   template<int nDim>
   std::ostream& operator<<( std::ostream& os, const Point<nDim>& p );

   template<int nDim>
   struct Direction
  {
      using Point = Geometry::Point<nDim>;

      std::array<Types::Real,nDim>  x{};

   // default, copy and move constructors
      Direction() = default;
      Direction( const Direction&  ) = default;
      Direction(       Direction&& ) = default;

   // only explicit construction from VariableDelta
      explicit Direction( const Point&  d ) noexcept : x(d.x) {}
      explicit Direction(       Point&& d ) noexcept : x(std::move(d.x)) {}

   // copy/move assignment
      Direction& operator=( const Direction&  ) = default;
      Direction& operator=(       Direction&& ) = default;

   // initialiser list constructor, must be a length nDim list of Types::Real
      template<typename... T>
         requires   Same<Types::Real,T...>
                 && is_pack_of_n<nDim,T...>::value
      Direction( T... r ) noexcept : x{r...} {}

   // accessors
            Types::Real& operator[]( const int i )       { return x[i]; };
      const Types::Real& operator[]( const int i ) const { return x[i]; }

   // in-place arithmetic
      Direction& operator+=( const Direction&  d );
      Direction& operator-=( const Direction&  d );
      Direction& operator*=( const Types::Real a );
      Direction& operator/=( const Types::Real a );
      Direction& operator =( const Types::Real a );
  };

/*
 * send each element of Direction to stream, seperated by a single space
 */
   template<int nDim>
   std::ostream& operator<<( std::ostream& os, const Direction<nDim>& d );

/*
 * Metric for a surface in nDim dimensions. Direction m[0] is surface normal, m[1...] are surface tangents
 */
   template<int nDim>
   struct Metric
  {
      std::array<Direction<nDim>,nDim>  m{};

            Direction<nDim>& operator[](       int i )       { return m[i]; };
      const Direction<nDim>& operator[]( const int i ) const { return m[i]; };
  };

// used to represent hypersurfaces in a space
   template<int nDim>
   struct Surface
  {
      Types::Real     area{};
      Point<nDim>   centre{};
      Metric<nDim>  metric{};
  };

// used to represent sections of a space of the same # of dimensions as the space
   template<int nDim>
   struct Volume
  {
      Types::Real  volume{};
      Point<nDim>  centre{};
  };

   template<int nDim>
   Types::Real length2( const Direction<nDim>& );

   template<int nDim>
   Types::Real length(  const Direction<nDim>& );

   Direction<2> cross( const Direction<2>& );

   Direction<3> cross( const Direction<3>&,
                       const Direction<3>& );

   Surface<1> surface( const Point<1>& );

   Surface<2> surface( const Point<2>&,
                       const Point<2>& );

   Surface<3> surface( const Point<3>&,
                       const Point<3>&,
                       const Point<3>&,
                       const Point<3>& );

   Volume<1> volume( const Point<1>&, const Point<1>& );

   Volume<2> volume( const Point<2>&, const Point<2>&,
                     const Point<2>&, const Point<2>& );

   Volume<3> volume( const Point<3>&, const Point<3>&,
                     const Point<3>&, const Point<3>&,
                     const Point<3>&, const Point<3>& );

   std::vector<Volume<1>> dual( const std::vector<Point<1>>&  nodes );
   std::vector<Point<1>>  dual( const std::vector<Volume<1>>& cells );


# include <geometry/direction.ipp>
# include <geometry/point.ipp>
# include <geometry/affine.ipp>
# include <geometry/operations.ipp>
# include <geometry/surface.ipp>
# include <geometry/volume.ipp>
# include <geometry/dual.ipp>

}

