
# include <vector>

# include <cassert>

namespace geom
{

   template<floating_point Real>
   std::vector<Volume<1,Real>> dual( const std::vector<Point<1,Real>>& nodes )
  {
      std::vector<Volume<1,Real>> cells(nodes.size()-1);
      dual( nodes, cells );
      return cells;
  }

   template<floating_point Real>
   void dual( const std::vector<Point< 1,Real>>& nodes,
                    std::vector<Volume<1,Real>>& cells )
  {
      const auto n=cells.size();
      assert( n+1 == nodes.size() );

      for( size_t i=0; i<n; i++ )
     {
         cells[i] = volume( nodes[i], nodes[i+1] );
     }
      return;
  }
}
