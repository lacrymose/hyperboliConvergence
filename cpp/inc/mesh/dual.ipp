
/*
 * return value form
 */
   template<int nDim, floating_point Real>
   MeshCellArray<nDim,Real> dual( const MeshNodeArray<nDim,Real>& nodes )
  {
      MeshCellArray<nDim,Real> cells(dualShape(nodes.shape()));
      dual( nodes, cells );
      return cells;
  }

/*
 * 1D dual calculation
 */
   template<floating_point Real>
   void dual( const MeshNodeArray<1,Real>& nodes,
                    MeshCellArray<1,Real>& cells )
  {
      assert( cells.shape() == par::dualShape( nodes.shape() ) );

      using Idx_n = typename MeshNodeArray<1,Real>::IdxType;
      using Idx_c = typename MeshCellArray<1,Real>::IdxType;

      const auto n = cells.shape(0);

      for( size_t i=0; i<n; ++i )
     {
         const Idx_n i_n0{i};
         const Idx_n i_n1{i+1};
         const Idx_c i_c{i};

         cells[i_c] = geom::volume( nodes[i_n0],
                                    nodes[i_n1] );
     }
      return;
  }

/*
 * 2D dual calculation
 */
   template<floating_point Real>
   void dual( const MeshNodeArray<2,Real>& nodes,
                    MeshCellArray<2,Real>& cells )
  {
      assert( cells.shape() == par::dualShape( nodes.shape() ) );

      using Idx_n = typename MeshNodeArray<2,Real>::IdxType;
      using Idx_c = typename MeshCellArray<2,Real>::IdxType;

      using Off_n = typename MeshNodeArray<2,Real>::OffsetType;

   // corners of standard quad
   //
   //   2 ---- 3
   //   |      |
   //   |      |
   //   0 ---- 1

      constexpr Off_n o0{0,0};
      constexpr Off_n o1{1,0};
      constexpr Off_n o2{0,1};
      constexpr Off_n o3{1,1};

      const auto ncx = cells.shape(0);
      const auto ncy = cells.shape(1);

      for( size_t i=0; i<ncx; ++i )
     {
         for( size_t j=0; j<ncy; ++j )
        {
            const Idx_c ij_c{i,j};
            const Idx_n ij_n{i,j};

            cells[ij_c] = geom::volume( nodes[ ij_n+o0 ],
                                        nodes[ ij_n+o1 ],
                                        nodes[ ij_n+o2 ],
                                        nodes[ ij_n+o3 ] );
        }
     }

      return;
  }

