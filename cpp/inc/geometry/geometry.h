
# pragma once

# include <geometry/declarations.h>
# include <geometry/type-traits.h>

# include <utils/maths/affineSpace.h>

# include <utils/concepts.h>

# include <array>

namespace geom
{
// ----------------- geometry structs ----------------- 

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


/*
 * used to represent hyper-quadrilaterals in an nDim space
 */
   template<int nDim, floating_point Real>
   struct Surface
  {
      Real                area;
      Point<nDim,Real>  centre;
      Metric<nDim,Real> metric;
  };


/*
 * used to represent hyper-hexahedrons in an nDim space
 */
   template<int nDim, floating_point Real>
   struct Volume
  {
      Real             volume;
      Point<nDim,Real> centre;
  };


// ----------------- functions on metrics ----------------- 

/*
 *  return identity metric
 */
   template<int nDim, floating_point Real>
   Metric<nDim,Real> metricI();

/*
 * return transpose metric
 *    metrics are unary matrices, so transpose is inverse
 */
   template<int nDim, floating_point Real>
   Metric<nDim,Real> transpose( const Metric<nDim,Real>& m );


// ----------------- operations on Points/Directions ----------------- 

/*
 * Dot product of two Directions
 */
   template<int nDim, floating_point Real>
   Real dot( const Direction<nDim,Real>&,
             const Direction<nDim,Real>& );

/*
 * Length of a Direction squared
 *    comparison of Direction lengths can be done on length2 with the same result, avoiding a sqrt
 */
   template<int nDim, floating_point Real>
   Real length2( const Direction<nDim,Real>& );

/*
 * Length of a Direction
 */
   template<int nDim, floating_point Real>
   Real length(  const Direction<nDim,Real>& );


/*
 * Cross products of two directions
 */
   template<floating_point Real>
   Real cross( const Direction<2,Real>&,
               const Direction<2,Real>& );

   template<floating_point Real>
   Direction<3,Real> cross( const Direction<3,Real>&,
                            const Direction<3,Real>& );

/*
 * Orthogonal vector in 2D
 */
   template<floating_point Real>
   Direction<2,Real> orthog( const Direction<2,Real>& );

/*
 * coefficients of equation describing straight line y=ax+b
 */
   template<floating_point Real>
   std::array<Real,2> line_coefficients( const Point<2,Real>& p, const Direction<2,Real>& d );

/*
 * vertex centroid of a polygon is average of all vertices
 */
   template<typename Pt0, typename... Pts>
      requires    is_Point_v<Pt0>
              && ((std::is_same_v<Pt0,Pts>)&&...)
   Pt0 vertex_centroid( const Pt0& pt0,  const Pts&... pts );

/*
 * mass centroid of a polygon is centre of mass
 */
   template<floating_point Real>
   Point<2,Real> mass_centroid( const Point<2,Real>& p0, const Point<2,Real>& p1,
                                const Point<2,Real>& p2, const Point<2,Real>& p3 );


// ----------------- generation of geometric entities from Points ----------------- 

/*
 * Create a surface (one less dimension than the background space) from its corners
 */
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


/*
 * Create a volume (same number of dimensions as the background space) from its corners
 */
   template<floating_point Real>
   Volume<1,Real> volume( const Point<1,Real>&, const Point<1,Real>& );

   template<floating_point Real>
   Volume<2,Real> volume( const Point<2,Real>&, const Point<2,Real>&,
                          const Point<2,Real>&, const Point<2,Real>& );

   template<floating_point Real>
   Volume<3,Real> volume( const Point<3,Real>&, const Point<3,Real>&,
                          const Point<3,Real>&, const Point<3,Real>&,
                          const Point<3,Real>&, const Point<3,Real>& );


}

# include <geometry/operations.ipp>

# include <geometry/metric.ipp>
# include <geometry/surface.ipp>
# include <geometry/volume.ipp>

# include <geometry/centroid.ipp>

