
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

   // direction from interior cell to boundary face
      const geom::Direction<nDim,Real> toBoundary = face.centre - cellinterior.centre;

   // face normal points into or out of domain?
      const bool face_into_domain = dot( toBoundary, face.metric[0] ) < 0;

   // face normal points to right value
      const SolVarT& qr = face_into_domain ? qinterior : qboundary;
      const SolVarT& ql = face_into_domain ? qboundary : qinterior;

   // face aligned average
      const StateT stateA = set2State( species,
                                       rotateToFace( face,
                                                     ql+0.5*(qr-ql) ) );

      if( stateA.velocity(0) > 0 ) // flow in direction of face normal
     {
         return ql;
     }
      else if( stateA.velocity(0) < 0 ) // flow in opposite direction to face normal
     {
         return ql;
     }
      else
     {
         return rotateFromFace( face,
                                state2Set<SolVarT>( species,
                                                    stateA ) );
     }
  }
