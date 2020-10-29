
# pragma once

# include <conservationLaws/euler/euler.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>

# include <type_traits>

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
                           const SolVarT&                    qboundary )
  {
      return qboundary;
  }
   
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
                     const SolVarT&                    qfarfield )
  {
      using FluxRes = FluxResult<LawType::Euler,nDim,Real>;
      using StateT  = State<LawType::Euler,nDim,Real>;

//    const int direction = (celli.centre[1]>0.5) ? 1 : -1;
//    const int direction = 0;
      const int direction = -1;
      const SolVarT qwall = qinterior + direction*dq;
      const StateT stateWall = set2State( species, qwall );
      const Real pressureWall = stateWall.pressure();

   // convective flux is zero (no flow through wall)
   // only flux is pressure flux
      FluxRes fr;
      for( unsigned int i=0; i<nDim; ++i )
     {
         fr.flux[i] = face.area*pressureWall*face.metric[0][i];
     }
      fr.flux[nDim  ] = 0.;
      fr.flux[nDim+1] = 0.;

      fr.lambda = face.area*spectralRadius( face.metric[0], stateWall );

      return fr;
  }
