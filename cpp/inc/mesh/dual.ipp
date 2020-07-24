
/*
 * return value form
 */
   template<int nDim, floating_point Real>
   par::Array<geom::Volume<nDim,Real>,nDim> dual( const par::Array<geom::Point<nDim,Real>,nDim>& nodes )
  {
      par::Array<geom::Volume<nDim,Real>,nDim> cells(cellDims_from_nodeDims(nodes.shape()));
      dual( nodes, cells );
      return cells;
  }

/*
 * 1D dual calculation
 */
   template<floating_point Real>
   void dual( const par::Array<geom::Point< 1,Real>,1>& nodes,
                    par::Array<geom::Volume<1,Real>,1>& cells )
  {
      const auto n=cells.shape(0);
      assert( n+1 == nodes.shape(0) );

      for( size_t i=0; i<n; i++ )
     {
         cells[{i}] = geom::volume( nodes[{i}], nodes[{i+1}] );
     }
      return;
  }

/*
 * 2D dual calculation
 */
   template<floating_point Real>
   void dual( const par::Array<geom::Point< 2,Real>,2>& nodes,
                    par::Array<geom::Volume<2,Real>,2>& cells )
  {
      assert( cells.shape() == cellDims_from_nodeDims(nodes.shape()) );

      const size_t ncx = cells.shape(0);
      const size_t ncy = cells.shape(1);

   // corners of standard quad
   //
   //   2 ---- 3
   //   |      |
   //   |      |
   //   0 ---- 1

      const par::Offset<2> o0{0,0};
      const par::Offset<2> o1{1,0};
      const par::Offset<2> o2{0,1};
      const par::Offset<2> o3{1,1};

      for( size_t i=0; i<ncx; i++ )
     {
         for( size_t j=0; j<ncy; j++ )
        {
            const par::Idx<2> ij{i,j};

            cells[ij] = geom::volume( nodes[ ij+o0 ], nodes[ ij+o1 ],
                                      nodes[ ij+o2 ], nodes[ ij+o3 ] );
        }
     }

      return;
  }

