
# pragma once

# include <conservationLaws/scalarAdvection/scalarAdvection.h>

/*
 * calculate boundary value from Riemann invariant conditions
 */
   template<ScalarVarSet SolVarT,
            int             nDim,
            floating_point  Real>
      requires SameDim<SolVarT,
                       dim_constant<nDim>>
   SolVarT boundaryUpdate(       std::integral_constant<ScalarAdvectionBCs,ScalarAdvectionBCs::Riemann>,
                           const Species<LawType::ScalarAdvection,Real>& species,
                           const geom::Surface<nDim,Real>&        face,
                           const geom::Volume<nDim,Real>&        cell,
                           const SolVarT&                           ql,
                           const SolVarT&                           qr );
   
# include <conservationLaws/scalarAdvection/boundary/riemann.ipp>

