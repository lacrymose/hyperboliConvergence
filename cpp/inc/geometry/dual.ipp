
# include <vector>

# include <cassert>

namespace geom
{

   template<floating_point Real>
   MDArray<Volume<1,Real>,1> dual( const MDArray<Point<1,Real>,1>& nodes )
  {
      MDArray<Volume<1,Real>,1> cells(nodes.dims[0]-1);
      dual( nodes, cells );
      return cells;
  }

   template<floating_point Real>
   void dual( const MDArray<Point< 1,Real>,1>& nodes,
                    MDArray<Volume<1,Real>,1>& cells )
  {
      const auto n=cells.dims[0];
      assert( n+1 == nodes.dims[0] );

      for( size_t i=0; i<n; i++ )
     {
         cells[{i}] = volume( nodes[{i}], nodes[{i+1}] );
     }
      return;
  }
}
