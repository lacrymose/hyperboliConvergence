
# pragma once

# include <geometry/geometry.h>

# include <parallalg/array.h>

   template<int            nDim,
            floating_point Real>
   using MeshNodeArray = par::PrimalArray<geom::Point<nDim,Real>,nDim>;

   template<int            nDim,
            floating_point Real>
   using MeshCellArray = par::DualArray<geom::Volume<nDim,Real>,nDim>;

/*
 * data structure holding the cell Volumes and node Points for a structured mesh
 */
   template<int            nDim,
            floating_point Real>
   struct Mesh
  {
      using Node = geom::Point<  nDim,Real>;
      using Cell = geom::Volume< nDim,Real>;
      using Face = geom::Surface<nDim,Real>;

      using NodeArray = par::PrimalArray<Node,nDim>;
      using CellArray = par::DualArray<  Cell,nDim>;

      using NodeShape = par::PrimalShape<nDim>;
      using CellShape = par::DualShape<  nDim>;

      NodeShape node_shape;
      CellShape cell_shape;

      NodeArray nodes;
      CellArray cells;

   // par::Array only supports move construction, so same must be for Mesh
      Mesh() = delete;
      Mesh( const Mesh&  ) = delete;
      Mesh(       Mesh&& ) = default;

   // par::Array only supports move assignment, so same must be for Mesh
      Mesh& operator=( const Mesh&  ) = delete;
      Mesh& operator=(       Mesh&& ) = default;

      Mesh( const NodeShape& s ) : node_shape(s),
                                   cell_shape(par::dualShape(s)),
                                   nodes(node_shape),
                                   cells(cell_shape) {}

      Mesh( const CellShape& s ) : node_shape(par::primalShape(s)),
                                   cell_shape(s),
                                   nodes(node_shape),
                                   cells(cell_shape) {}
  };

// ----------------- operations on the dual mesh ----------------- 

/*
 * Create the dual mesh from the primal mesh (assumes mesh is structured)
 *    the primal mesh is the nodes
 *    the dual mesh is the cells
 */
   template<int nDim, floating_point Real>
   MeshCellArray<nDim,Real> dual( const MeshNodeArray<nDim,Real>& nodes );

   template<floating_point Real>
   void dual( const MeshNodeArray<1,Real>& nodes,
                    MeshCellArray<1,Real>& cells );

   template<floating_point Real>
   void dual( const MeshNodeArray<2,Real>& nodes,
                    MeshCellArray<2,Real>& cells );


# include <mesh/dual.ipp>

