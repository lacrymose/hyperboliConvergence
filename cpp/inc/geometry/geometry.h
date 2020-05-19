
# pragma once

# include <utils/maths/affineSpace.h>

# include <utils/concepts.h>

# include <array>
# include <vector>

namespace Geometry
{
// forward declarations
   template<int nDim, floating_point Real=double> struct Point;
   template<int nDim, floating_point Real=double> struct Direction;
   template<int nDim, floating_point Real=double> struct Metric;
   template<int nDim, floating_point Real=double> struct Volume;
   template<int nDim, floating_point Real=double> struct Surface;

/*
 * {Point,Direction} form an affine space.
 *    affine space behaviour is inherited from CRTP base classes
 */
   template<int nDim, floating_point Real>
   struct Point : AffinePointBase<nDim,
                                  Point<nDim,Real>,
                                  Direction<nDim,Real>,
                                  Real>
  {
      using AffinePointBase<nDim,
                            Point<nDim,Real>,
                            Direction<nDim,Real>,
                            Real>::AffinePointBase;
  };

   template<int nDim, floating_point Real>
   struct Direction : AffineDeltaBase<nDim,
                                      Point<nDim,Real>,
                                      Direction<nDim,Real>,
                                      Real>
  {
      using AffineDeltaBase<nDim,
                            Point<nDim,Real>,
                            Direction<nDim,Real>,
                            Real>::AffineDeltaBase;
  };

/*
 * Metric for a surface in nDim dimensions. Direction m[0] is surface normal, m[1...] are surface tangents
 */
   template<int nDim, floating_point Real>
   struct Metric
  {
      std::array<Direction<nDim,Real>,nDim>  m;

            Direction<nDim,Real>& operator[](       int i )       { return m[i]; };
      const Direction<nDim,Real>& operator[]( const int i ) const { return m[i]; };
  };

// used to represent hypersurfaces in a space
   template<int nDim, floating_point Real>
   struct Surface
  {
      Real                area;
      Point<nDim,Real>  centre;
      Metric<nDim,Real> metric;
  };

// used to represent sections of a space of the same # of dimensions as the space
   template<int nDim, floating_point Real>
   struct Volume
  {
      Real             volume;
      Point<nDim,Real> centre;
  };

   template<int nDim, floating_point Real>
   Real length2( const Direction<nDim,Real>& );

   template<int nDim, floating_point Real>
   Real length(  const Direction<nDim,Real>& );

   template<floating_point Real>
   Direction<2,Real> cross( const Direction<2,Real>& );

   template<floating_point Real>
   Direction<3,Real> cross( const Direction<3,Real>&,
                            const Direction<3,Real>& );

   template<floating_point Real>
   Surface<1,Real> surface( const Point<1,Real>& );

   template<floating_point Real>
   Surface<2,Real> surface( const Point<2,Real>&,
                            const Point<2,Real>& );

   template<floating_point Real>
   Surface<3,Real> surface( const Point<3,Real>&,
                            const Point<3,Real>&,
                            const Point<3,Real>&,
                            const Point<3,Real>& );

   template<floating_point Real>
   Volume<1,Real> volume( const Point<1,Real>&, const Point<1,Real>& );

   template<floating_point Real>
   Volume<2,Real> volume( const Point<2,Real>&, const Point<2,Real>&,
                          const Point<2,Real>&, const Point<2,Real>& );

   template<floating_point Real>
   Volume<3,Real> volume( const Point<3,Real>&, const Point<3,Real>&,
                          const Point<3,Real>&, const Point<3,Real>&,
                          const Point<3,Real>&, const Point<3,Real>& );

   template<floating_point Real>
   std::vector<Volume<1,Real>> dual( const std::vector<Point<1,Real>>& nodes );

   template<floating_point Real>
   std::vector<Point<1,Real>> dual( const std::vector<Volume<1,Real>>& cells );
}

# include <geometry/operations.ipp>
# include <geometry/surface.ipp>
# include <geometry/volume.ipp>
# include <geometry/dual.ipp>

