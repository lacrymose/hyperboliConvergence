
# pragma once

# include <mesh/mesh.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <array>
# include <cmath>

/*
 * create a 2D mesh over domain [lo0:hi0]x[lo1:hi1] with linear spacing
 */
   template<floating_point Real>
   Mesh<2,Real> make_linspace_mesh( const par::Shape<2>& cellDims, const Real lox, const Real hix,
                                                                   const Real loy, const Real hiy )
  {
      Mesh<2,Real> mesh(cellDims);

   // number of cells in each direction
      const size_t ncx = mesh.cells.shape(0);
      const size_t ncy = mesh.cells.shape(1);

      const Real dx = ( hix - lox ) / ncx;
      const Real dy = ( hiy - loy ) / ncy;

      using Pt = typename Mesh<2,Real>::Node;
      par::generate_idx( mesh.nodes,
                         // return linear spacing in x and y
                          [=]( par::Idx<2> i ) -> Pt
                         { return Pt{lox+i[0]*dx,loy+i[1]*dy}; }
                       );

   // initialise cell array
      mesh.cells = dual( mesh.nodes );

      return mesh;
  }

/*
 * create a 2D mesh of a channel with a sinusoidal bump
 *    d is height of channel
 *    l is length of channel before/after bump
 *    w is width of bump
 *    h is height of bump
 */
   template<floating_point Real>
   Mesh<2,Real> make_channelbump_mesh( const std::array<Real,4>  lengths,
                                       const std::array<size_t,3> ncells )
  {
      using Node = typename Mesh<2,Real>::Node;

      const Real d = lengths[0];
      const Real l = lengths[1];
      const Real w = lengths[2];
      const Real h = lengths[3];

   // number of cells
      const size_t nd = ncells[0];
      const size_t nl = ncells[1];
      const size_t nw = ncells[2];

      const size_t nx = 2*nl+nw;
      const size_t ny = nd;

      Mesh<2,Real> mesh( par::Shape<2>{nx,ny} );
      assert( mesh.nodes.shape(0) == 2*nl+nw+1 );
      assert( mesh.nodes.shape(1) ==      ny+1 );

   // x coordinate at left edge of left/right channel region and bump
      const Real xl0 = -( l + 0.5*w );
      const Real xw0 =       -0.5*w;
      const Real xr0 =        0.5*w;

   // nodes before and after bump
      const Real dxl = l/nl;
      const Real dyl = d/ny;
      for( size_t i=0; i<nl+1; ++i )
     {
      // offset for downstream of bump
         const size_t off = nl+nw;

         for( size_t j=0; j<ny+1; ++j )
        {
            const Real xin  = xl0 + i*dxl;
            const Real xout = xr0 + i*dxl;
            const Real y = j*dyl;

            const Node pt0{ xin , y };
            const Node pt1{ xout, y };

            mesh.nodes[{     i, j }] = pt0;
            mesh.nodes[{ off+i, j }] = pt1;
        }
     }

   // cells over bump
      const Real dxw = w/nw;
      for( size_t i=0; i<nw; ++i )
     {
         const size_t ii = i+nl+1;

         const Real x = xw0 + (i+1)*dxw;

      // node on lower (bump) and higher (flat) walls
         const Real yl = 0.5*h*( 1. + cos( 2.*M_PI*x/w ) );
         const Real yh = d;

         const Node pl = Node{ x, yl };
         const Node ph = Node{ x, yh };

         mesh.nodes[{ ii, 0  }] = pl;
         mesh.nodes[{ ii, ny }] = ph;

      // delta between lower and higher walls
         const auto dp = (ph-pl)*(1./ny);

      // linear distribution from bottom to top
         for( size_t j=0; j<ny; ++j )
        {
            mesh.nodes[{ ii, j }] = pl + Real(j)*dp;
        }
     }

   // initialise cell array
      mesh.cells = dual( mesh.nodes );

      return mesh;
  }

