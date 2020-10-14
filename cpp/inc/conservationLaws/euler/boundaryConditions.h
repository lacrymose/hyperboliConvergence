
# pragma once

# include <conservationLaws/euler/euler.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>
# include <type_traits>

// ---------------------- Riemann invariant ----------------------

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
                           const geom::Volume<nDim,Real>&        celli,
                           const SolVarT&                    qinterior,
                           const SolVarT&                    qboundary );
   
// ---------------------- Entropy freestream ----------------------

/*
 * calculate boundary value, assigning velocity and entropy at inlet and pressure at outlet
 */
   template<EulerVarSet SolVarT,
            int            nDim,
            floating_point Real>
      requires SameDim<SolVarT,
                       dim_constant<nDim>>
   SolVarT boundaryUpdate(       std::integral_constant<EulerBCs,EulerBCs::Entropy>,
                           const Species<LawType::Euler,Real>& species,
                           const geom::Surface<nDim,Real>&        face,
                           const geom::Volume<nDim,Real>& cellinterior,
                           const SolVarT&                    qinterior,
                           const SolVarT&                    qfarfield );

// ---------------------- Inviscid wall ----------------------

/*
 * calculate boundary value for an inviscid wall (boundary value unused in flux calculation)
 */
   template<EulerVarSet SolVarT,
            int            nDim,
            floating_point Real>
      requires SameDim<SolVarT,
                       dim_constant<nDim>>
   SolVarT boundaryUpdate(       std::integral_constant<EulerBCs,EulerBCs::InviscidWall>,
                           const Species<LawType::Euler,Real>& species,
                           const geom::Surface<nDim,Real>&        face,
                           const geom::Volume<nDim,Real>&        celli,
                           const SolVarT&                    qinterior,
                           const SolVarT&                    qboundary );
   
/*
 * calculate boundary flux for an inviscid wall
 */
   template<EulerVarSet   SolVarT,
            EulerVarDelta SolDelT,
            int              nDim,
            floating_point   Real>
      requires SameDim<SolVarT,
                       dim_constant<nDim>>
   FluxResult<LawType::Euler,nDim,Real>
       boundaryFlux(       std::integral_constant<EulerBCs,EulerBCs::InviscidWall>,
                     const Species<LawType::Euler,Real>& species,
                     const geom::Surface<nDim,Real>&        face,
                     const geom::Volume<nDim,Real>&        celli,
                     const SolDelT                            dq,
                     const SolVarT&                    qinterior,
                     const SolVarT&                    qboundary,
                     const SolVarT&                    qfarfield );

# include <conservationLaws/euler/boundary/riemann.ipp>
# include <conservationLaws/euler/boundary/entropy.ipp>
# include <conservationLaws/euler/boundary/inviscidWall.ipp>

