
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>

# include <numeric>
# include <algorithm>
# include <limits>

# include <cassert>

   template<floating_point Real, int nDim, LawType Law>
   Real spectralRadius( const std::vector<geom::Volume<  nDim,Real>>& cells,
                        const std::vector<FluxResult<Law,nDim,Real>>& resid )
  {
      assert( cells.size() == resid.size() );

   // max( (cell,resid)->spectral radius ... )
      return std::inner_product( cells.begin(), cells.end(),
                                 resid.begin(),

                              // initial value for finding maximum
                                 std::numeric_limits<Real>::min(),

                              // find maximum of all spectral radii
                                 []( const Real l, const Real r )
                                { return std::max( l,r ); }

                              // find spectral radii of one cell/resid pair
                                 []( const geom::Volume<  nDim,Real>& cell,
                                     const FluxResult<Law,nDim,Real>& res )
                                { return res.lambda / cell.volume; }

                               );
  }

