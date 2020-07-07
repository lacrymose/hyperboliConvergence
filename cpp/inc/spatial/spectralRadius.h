
# pragma once

# include <conservationLaws/base/base.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <numeric>
# include <algorithm>
# include <limits>

# include <cassert>

   template<floating_point Real, int nDim, LawType Law>
   Real spectralRadius( const par::Array<geom::Volume<  nDim,Real>,1>& cells,
                        const par::Array<FluxResult<Law,nDim,Real>,1>& resid )
  {
      assert( cells.shape() == resid.shape() );

//    return std::inner_product( cells.elems.begin(), cells.elems.end(),
//                               resid.elems.begin(),

//                            // initial value for finding maximum
//                               std::numeric_limits<Real>::min(),

//                            // find maximum of all spectral radii
//                               []( const Real l, const Real r )
//                              { return std::max( l,r ); },

//                            // find spectral radii of one cell/resid pair
//                               []( const geom::Volume<  nDim,Real>& cell,
//                                   const FluxResult<Law,nDim,Real>& res )
//                              { return res.lambda / cell.volume; }

//                             );

   // max( (cell,resid)->spectral radius ... )
      return par::transform_reduce(
                                 // transform one cell/resid pair to a spectral radius
                                    []( const geom::Volume<  nDim,Real>& cell,
                                        const FluxResult<Law,nDim,Real>& res )
                                   { return res.lambda / cell.volume; },

                                 // reduce to maximum of all spectral radii
                                    []( const Real l, const Real r )
                                   { return std::max( l,r ); },

                                 // initial value for finding maximum
                                    std::numeric_limits<Real>::min(),

                                 // source arrays
                                    cells,
                                    resid
                                  );

  }

