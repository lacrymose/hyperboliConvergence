
   std::vector<Volume<1>> dual( const std::vector<Point<1>>& nodes )
  {
      const auto n=nodes.size();
      std::vector<Volume<1>> cells(n-1);

      for( size_t i=0; i<n-1; i++ )
     {
         cells[i] = volume( nodes[i], nodes[i+1] );
     }
      return cells;
  }
