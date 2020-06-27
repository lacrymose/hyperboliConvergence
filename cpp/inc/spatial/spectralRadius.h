
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>

# include <mdarray/mdarray.h>

# include <numeric>
# include <algorithm>
# include <limits>

# include <cassert>

   template<floating_point Real, int nDim, LawType Law>
   Real spectralRadius( const MDArray<geom::Volume<  nDim,Real>,1>& cells,
                        const MDArray<FluxResult<Law,nDim,Real>,1>& resid )
  {
      assert( cells.dims == resid.dims );

   // max( (cell,resid)->spectral radius ... )
      return std::inner_product( cells.elems.begin(), cells.elems.end(),
                                 resid.elems.begin(),

                              // initial value for finding maximum
                                 std::numeric_limits<Real>::min(),

                              // find maximum of all spectral radii
                                 []( const Real l, const Real r )
                                { return std::max( l,r ); },

                              // find spectral radii of one cell/resid pair
                                 []( const geom::Volume<  nDim,Real>& cell,
                                     const FluxResult<Law,nDim,Real>& res )
                                { return res.lambda / cell.volume; }

                               );
  }

