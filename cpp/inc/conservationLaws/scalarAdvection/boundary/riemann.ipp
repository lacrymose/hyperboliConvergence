
# pragma once

# include <conservationLaws/scalarAdvection/scalarAdvection.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>

# include <type_traits>

/*
 * calculate boundary value from Riemann invariant conditions
 */
   template<ScalarVarSet SolVarT,
            int             nDim,
            floating_point  Real>
      requires SameDim<SolVarT,
                       dim_constant<nDim>>
   SolVarT boundaryUpdate(       std::integral_constant<ScalarAdvectionBCs,
                                                        ScalarAdvectionBCs::Riemann>,
                           const Species<LawType::ScalarAdvection,Real>& species,
                           const geom::Surface<nDim,Real>&                  face,
                           const geom::Volume<nDim,Real>&           cellinterior,
                           const SolVarT&                              qinterior,
                           const SolVarT&                              qboundary )
  {
      using StateT = State<LawType::ScalarAdvection,nDim,Real>;

   // face aligned average
      const StateT stateA = set2State( species,
                                       rotateToFace( face,
                                                     qinterior+0.5*(qboundary-qinterior) ) );

   // face normal points into domain
      const bool inflow = projectedVelocity( face, stateA )>0;

      if    (  inflow )  { return qboundary; }
      else /* outflow */ { return qinterior; }
  }
