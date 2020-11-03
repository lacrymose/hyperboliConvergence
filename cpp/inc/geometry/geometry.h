
# pragma once

# include <geometry/declarations.h>
# include <geometry/type-traits.h>

# include <utils/maths/affineSpace.h>

# include <utils/concepts.h>

# include <array>

namespace geom
{
/*
 * ----------------- geometry structs ----------------- 
 */

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

/*
 * Global Point object at coordinate origin
 */
   template<int nDim, floating_point Real>
   constexpr Point<nDim,Real> origin{};

   template<floating_point Real>
   constexpr Point<1,Real> origin1{};

   template<floating_point Real>
   constexpr Point<2,Real> origin2{};

   template<floating_point Real>
   constexpr Point<3,Real> origin3{};

/*
 * ----------------- Convenience typedefs for geometry structs ------------
 */

// Point
   template<floating_point Real>
   using Point1 = Point<1,Real>;

   template<floating_point Real>
   using Point2 = Point<2,Real>;

   template<floating_point Real>
   using Point3 = Point<3,Real>;

// Direction
   template<floating_point Real>
   using Direction1 = Direction<1,Real>;

   template<floating_point Real>
   using Direction2 = Direction<2,Real>;

   template<floating_point Real>
   using Direction3 = Direction<3,Real>;

// Metric
   template<floating_point Real>
   using Metric1 = Metric<1,Real>;

   template<floating_point Real>
   using Metric2 = Metric<2,Real>;

   template<floating_point Real>
   using Metric3 = Metric<3,Real>;

// Surface
   template<floating_point Real>
   using Surface1 = Surface<1,Real>;

   template<floating_point Real>
   using Surface2 = Surface<2,Real>;

   template<floating_point Real>
   using Surface3 = Surface<3,Real>;

// Volume
   template<floating_point Real>
   using Volume1 = Volume<1,Real>;

   template<floating_point Real>
   using Volume2 = Volume<2,Real>;

   template<floating_point Real>
   using Volume3 = Volume<3,Real>;


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
 * Normalise direction to unit length
 */
   template<int nDim, floating_point Real>
   Direction<nDim,Real> norm( const Direction<nDim,Real>& d )

/*
 * Cross products of two directions
 */
   template<floating_point Real>
   Real cross( const Direction2<Real>&,
               const Direction2<Real>& );

   template<floating_point Real>
   Direction3<Real> cross( const Direction3<Real>&,
                           const Direction3<Real>& );

/*
 * Orthogonal vector in 2D
 */
   template<floating_point Real>
   Direction2<Real> orthog( const Direction2<Real>& );

/*
 * coefficients of equation describing straight line y=ax+b
 */
   template<floating_point Real>
   std::array<Real,2> line_coefficients( const Point2<Real>&     p,
                                         const Direction2<Real>& d );

/*
 * vertex centroid of a polygon is average of all vertices
 */
   template<typename Pt0, typename... Pts>
      requires    is_Point_v<Pt0>
              && ((std::is_same_v<Pt0,Pts>)&&...)
   Pt0 vertex_centroid( const Pt0&    pt0,
                        const Pts&... pts );

/*
 * mass centroid of a polygon is centre of mass
 */
   template<floating_point Real>
   Point2<Real> mass_centroid( const Point2<Real>& p0, const Point2<Real>& p1,
                               const Point2<Real>& p2, const Point2<Real>& p3 );

/*
 * Flip a direction to face the opposite way
 */
   template<int nDim, floating_point Real>
   Direction<nDim,Real> flip ( const Direction<nDim,Real>& d0 );

/*
 * Flip a surface so normal faces opposite way (maintaining right-handedness of coordinate metric)
 */
   template<int nDim, floating_point Real>
   Surface<nDim,Real> flip ( const Surface<nDim,Real>& s0 );



// ----------------- generation of geometric entities from Points ----------------- 

/*
 * Create a surface (one less dimension than the background space) from its corners
 */
   template<floating_point Real>
   Surface1<Real> surface( const Point1<Real>& );

   template<floating_point Real>
   Surface2<Real> surface( const Point2<Real>&,
                           const Point2<Real>& );

   template<floating_point Real>
   Surface3<Real> surface( const Point3<Real>&,
                           const Point3<Real>&,
                           const Point3<Real>&,
                           const Point3<Real>& );


/*
 * Create a volume (same number of dimensions as the background space) from its corners
 */
   template<floating_point Real>
   Volume1<Real> volume( const Point1<Real>&, const Point1<Real>& );
                                                          
   template<floating_point Real>                          
   Volume2<Real> volume( const Point2<Real>&, const Point2<Real>&,
                         const Point2<Real>&, const Point2<Real>& );
                                                          
   template<floating_point Real>                          
   Volume3<Real> volume( const Point3<Real>&, const Point3<Real>&,
                         const Point3<Real>&, const Point3<Real>&,
                         const Point3<Real>&, const Point3<Real>& );


}

# include <geometry/operations.ipp>

# include <geometry/metric.ipp>
# include <geometry/surface.ipp>
# include <geometry/volume.ipp>

# include <geometry/centroid.ipp>

