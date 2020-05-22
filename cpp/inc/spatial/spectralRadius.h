
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>

# include <numeric>
# include <algorithm>
# include <limits>

# include <cassert>

   template<floating_point Real, int nDim, LawType Law>
   Real spectralRadius( const std::vector<geom::Volume<nDim,Real>>& cells,
                        const std::vector<FluxResult<Law,nDim,Real>>&   resid )
  {
      assert( cells.size() == resid.size() );

   // functor to calculate spectral radius of single cell
      auto srCalc = []( const geom::Volume<nDim,Real>& cell,
                        const FluxResult<Law,nDim,Real>&   res )
     {
         return res.lambda / cell.volume;
     };

   // functor to accumulate the maximum spectral radius
      auto max = []( const Real l, const Real r )
     {
         return std::max( l,r );
     };

   // find max spectral radius
      return std::inner_product( cells.begin(), cells.end(),
                                 resid.begin(),
                                 std::numeric_limits<Real>::min(),
                                 max, srCalc );
  }

