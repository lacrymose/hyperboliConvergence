
# pragma once

# include <geometry/geometry.h>

# include <parallalg/array.h>

/*
 * data structure holding the cell Volumes and node Points for a structured mesh
 */
   template<int nDim, floating_point Real>
   struct Mesh
  {
      using Node = typename geom::Point<  nDim,Real>;
      using Cell = typename geom::Volume< nDim,Real>;
      using Face = typename geom::Surface<nDim,Real>;

      par::Shape<nDim> shape;

      par::Array<Node,nDim> nodes;
      par::Array<Cell,nDim> cells;

   // par::Array only supports move construction, so same must be for Mesh
      Mesh() = delete;
      Mesh( const Mesh&  ) = delete;
      Mesh(       Mesh&& ) = default;

   // par::Array only supports move assignment, so same must be for Mesh
      Mesh& operator=( const Mesh&  ) = delete;
      Mesh& operator=(       Mesh&& ) = default;

      Mesh( const par::Shape<nDim>& s ) : shape(s),
                                          nodes(par::nodeDims_from_cellDims(s)),
                                          cells(s) {}
  };

// ----------------- operations on the dual mesh ----------------- 

/*
 * Create the dual mesh from the primal mesh (assumes mesh is structured)
 *    the primal mesh is the nodes
 *    the dual mesh is the cells
 */
   template<int nDim, floating_point Real>
   par::Array<geom::Volume<nDim,Real>,nDim> dual( const par::Array<geom::Point<nDim,Real>,nDim>& nodes );

   template<floating_point Real>
   void dual( const par::Array<geom::Point< 1,Real>,1>& nodes,
                    par::Array<geom::Volume<1,Real>,1>& cells );

   template<floating_point Real>
   void dual( const par::Array<geom::Point< 2,Real>,2>& nodes,
                    par::Array<geom::Volume<2,Real>,2>& cells );


# include <mesh/dual.ipp>

