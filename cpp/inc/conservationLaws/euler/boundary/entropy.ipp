
# pragma once

# include <conservationLaws/euler/euler.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>

# include <array>
# include <type_traits>
# include <iostream>

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
                           const SolVarT&                    qfarfield )
  {
      using StateT = State<LawType::Euler,nDim,Real>;
      using PrimVarT = VariableSet<LawType::Euler,nDim,EulerBases::Primitive,Real>;
      using Direction = geom::Direction<nDim,Real>;

      const Real gam =species.gamma;

   // direction from boundary face into domain
      const Direction intoDomain = cellinterior.centre - face.centre;
   // face normal points into or out of domain?
      const bool face_into_domain = dot( intoDomain, face.metric[0] ) > 0;
      assert( face_into_domain );

   // interior/exterior/average states aligned with boundary face
      const StateT stateA = roeAverage( species,
                                        rotateToFace( face, qinterior ),
                                        rotateToFace( face, qfarfield ) );

   // if transverse flow, just return interior value
      if( ( stateA.velocity(0)
           / std::sqrt(stateA.velocity2()) ) < 0.0001 ){ return qinterior; }

   // boundary mach number
      const Real ma = stateA.velocity(0) / sqrt( stateA.speedOfSound2() );

      const bool inflow = ma > 0;

   // Upwind variables
      const SolVarT& quw = inflow ? qfarfield : qinterior;
      const SolVarT& qdw = inflow ? qinterior : qfarfield;

   // supersonic cases, return fully upwinded value
      if( std::abs(ma)>1. ){ return quw; }

   // subsonic cases:

   // interior/farfield states for linear fields
      const StateT stateU = set2State( species, rotateToFace( face, quw ) );
      const StateT stateD = set2State( species, rotateToFace( face, qdw ) );

      // pressure is downwind pressure
      const Real pb = stateD.pressure();

      // entropy is upwind entropy
      const Real sb =   stateU.speedOfSound2()
                     /( gam*pow( stateU.density(), gam-1. ) );

      // boundary density
      const Real rb = pow( pb/sb, 1./gam );

      PrimVarT qb;

      // velocity is upwind velocity
      for( int i=0; i<nDim; ++i ){ qb[i]=stateU.velocity(i); }
      qb[nDim]   = rb;
      qb[nDim+1] = pb;

//    return qfarfield;
      return rotateFromFace( face, set2Set<SolVarT>( species, qb ) );
  }
