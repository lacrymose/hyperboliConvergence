
# pragma once

# include <mesh/mesh.h>

# include <parallalg/array.h>

/*
 * create a 1D mesh over domain [lo:hi] with linear spacing and nc cells
 */
   template<floating_point Real>
   Mesh<1,Real> make_linspace_mesh( const par::DualShape1& cellDims, const Real lo, const Real hi )
  {
      Mesh<1,Real> mesh(cellDims);

      const size_t nc = cellDims[0];
      const size_t np = cellDims[0]+1;

      const Real dx = ( hi - lo ) / nc;

      using Pt = typename Mesh<1,Real>::Node;

      for( size_t i=0; i<np; i++ ){ mesh.nodes({i}) = Pt{lo+i*dx}; }

/*
      par::generate_idx( mesh.nodes,
                         // return linear spacing in x and y
                         [=]( par::PrimalIdx1 i ) -> Pt
                        { return Pt{lo+i[0]*dx}; }
                       );
*/

      mesh.cells = dual( mesh.nodes );

      return mesh;
  }

