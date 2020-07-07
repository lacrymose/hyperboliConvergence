
# pragma once

# include <utils/maths/affineSpace.h>

# include <utils/concepts.h>

# include <parallalg/array.h>

# include <array>

namespace geom
{
// forward declarations
   template<int nDim, floating_point Real> struct Point;
   template<int nDim, floating_point Real> struct Direction;
   template<int nDim, floating_point Real> struct Metric;
   template<int nDim, floating_point Real> struct Volume;
   template<int nDim, floating_point Real> struct Surface;

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
 * Struct holding the cell Volumes and node Points for a structured mesh
 */
   template<int nDim, floating_point Real>
   struct Mesh
  {
      par::Shape<nDim> shape;

      par::Array<Point< nDim,Real>,nDim> nodes;
      par::Array<Volume<nDim,Real>,nDim> cells;

   // par::Array only supports move construction, so same must be for Mesh
      Mesh() = delete;
      Mesh( const Mesh&  ) = delete;
      Mesh(       Mesh&& ) =default;

   // par::Array only supports move assignment, so same must be for Mesh
      Mesh& operator=( const Mesh&  ) = delete;
      Mesh& operator=(       Mesh&& ) = default;

      Mesh( const par::Shape<nDim>& s ) : shape(s),
                                     nodes(par::nodeDims_from_cellDims(s)),
                                     cells(s) {}
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
   par::Array<Volume<1,Real>,1> dual( const par::Array<Point<1,Real>,1>& nodes );

   template<floating_point Real>
   void dual( const par::Array<Point< 1,Real>,1>& nodes,
                    par::Array<Volume<1,Real>,1>& cells );

/*
 * create a 1D mesh over domain [lo:hi] with linear spacing and nc cells
 */
   template<floating_point Real>
   Mesh<1,Real> make_linspace_mesh( const par::Shape<1>& cellDims, const Real lo, const Real hi )
  {
      Mesh<1,Real> mesh(cellDims);

      const size_t nc = cellDims[0];
      const size_t np = cellDims[0]+1;

      const Real dx = ( hi - lo ) / nc;

      for( size_t i=0; i<np; i++ ){ mesh.nodes[{i}] = lo + i*dx; }
      mesh.cells = dual( mesh.nodes );

      return mesh;
  }
}

# include <geometry/operations.ipp>
# include <geometry/surface.ipp>
# include <geometry/volume.ipp>
# include <geometry/dual.ipp>

