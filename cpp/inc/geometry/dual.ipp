
# include <vector>

# include <cassert>

namespace geom
{

   template<floating_point Real>
   par::Array<Volume<1,Real>,1> dual( const par::Array<Point<1,Real>,1>& nodes )
  {
      par::Array<Volume<1,Real>,1> cells(nodes.shape(0)-1);
      dual( nodes, cells );
      return cells;
  }

   template<floating_point Real>
   void dual( const par::Array<Point< 1,Real>,1>& nodes,
                    par::Array<Volume<1,Real>,1>& cells )
  {
      const auto n=cells.shape(0);
      assert( n+1 == nodes.shape(0) );

      for( size_t i=0; i<n; i++ )
     {
         cells[{i}] = volume( nodes[{i}], nodes[{i+1}] );
     }
      return;
  }
}
