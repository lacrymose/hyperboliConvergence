
# pragma once

# include <mesh/mesh.h>

# include <parallalg/array.h>

/*
 * create a 2D mesh over domain [lo0:hi0]x[lo1:hi1] with linear spacing
 */
   template<floating_point Real>
   Mesh<2,Real> make_linspace_mesh( const par::Shape<2>& cellDims, const Real lox, const Real hix,
                                                                   const Real loy, const Real hiy )
  {
      Mesh<2,Real> mesh(cellDims);

   // number of nodes in each direction
      const size_t npx = cellDims[0]+1;
      const size_t npy = cellDims[0]+1;

   // number of cells in each direction
      const size_t ncx = cellDims[0];
      const size_t ncy = cellDims[0];

      const Real dx = ( hix - lox ) / ncx;
      const Real dy = ( hiy - loy ) / ncy;

      for( size_t i=0; i<npx; i++ )
     {
         for( size_t j=0; j<npy; j++ )
        {
            mesh.nodes[{i,j}][0] = lox + i*dx;
            mesh.nodes[{i,j}][1] = loy + j*dy;
        }
     }
      mesh.cells = dual( mesh.nodes );

/*
      using Pt = Point<2,Real>;
      par::generate_by_index( [=]( par::Idx<2> i )
                             {
                                 return Pt{lox+i[0]*dx,loy+i[1]*dy};
                             },
                              mesh.nodes
                            );
      mesh.cells = dual( mesh.nodes );
*/

      return mesh;
  }

