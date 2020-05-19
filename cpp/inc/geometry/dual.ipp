
# include <vector>

namespace Geometry
{
   template<floating_point Real>
   std::vector<Volume<1,Real>> dual( const std::vector<Point<1,Real>>& nodes )
  {
      const auto n=nodes.size();
      std::vector<Volume<1,Real>> cells(n-1);

      for( size_t i=0; i<n-1; i++ )
     {
         cells[i] = volume( nodes[i], nodes[i+1] );
     }
      return cells;
  }
}
