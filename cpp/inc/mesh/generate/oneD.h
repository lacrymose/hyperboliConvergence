
# pragma once

# include <mesh/mesh.h>

# include <parallalg/array.h>

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

      for( size_t i=0; i<np; i++ ){ mesh.nodes[{i}][0] = lo + i*dx; }
      mesh.cells = dual( mesh.nodes );

/*
      par::generate_by_index( [=]( par::Idx<1> i ){ return lo+i[0]*dx; },
                              mesh.nodes
                            );
      mesh.cells = dual( mesh.nodes );
*/

      return mesh;
  }

