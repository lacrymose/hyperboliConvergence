
# pragma once

# include <conservationLaws/euler/euler.h>

/*
 * calculate boundary value from Riemann invariant conditions
 */
   template<EulerVarSet SolVarT,
            int            nDim,
            floating_point Real>
      requires SameDim<SolVarT,
                       dim_constant<nDim>>
   SolVarT boundaryUpdate(       std::integral_constant<EulerBCs,EulerBCs::Riemann>,
                           const Species<LawType::Euler,Real>& species,
                           const geom::Surface<nDim,Real>&        face,
                           const SolVarT&                           ql,
                           const SolVarT&                           qr );
   
# include <conservationLaws/euler/boundary/riemann.ipp>

