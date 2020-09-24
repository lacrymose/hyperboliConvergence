
# pragma once

# include <conservationLaws/euler/euler.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>

# include <array>
# include <type_traits>
# include <iostream>

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
                           const geom::Volume<nDim,Real>& cellinterior,
                           const SolVarT&                    qinterior,
                           const SolVarT&                    qfarfield )
  {
      using StateT = State<LawType::Euler,nDim,Real>;
      using PrimVarT = VariableSet<LawType::Euler,nDim,EulerBases::Primitive,Real>;
      using Direction = geom::Direction<nDim,Real>;

      const Real gam =species.gamma;
      const Real gam1=species.gamma1;

   // direction from boundary face into domain
      const Direction intoDomain = cellinterior.centre - face.centre;
   // face normal points into or out of domain?
      const bool face_into_domain = dot( intoDomain, face.metric[0] ) > 0;
      assert( face_into_domain );

   // interior/exterior/average states aligned with boundary face
      const StateT stateA = roeAverage( species,
                                        rotateToFace( face, qinterior ),
                                        rotateToFace( face, qfarfield ) );

   // boundary mach number
      const Real ma = stateA.velocity(0) / sqrt( stateA.speedOfSound2() );

      const bool inflow = ma > 0;

   // Upwind variables
      const SolVarT& quw = inflow ? qfarfield : qinterior;

   // supersonic cases, return fully upwinded value
      if( std::abs(ma)>1. ){ return quw; }

   // subsonic cases:

   // interior/farfield states for linear fields
      const StateT stateI = set2State( species, rotateToFace( face, qinterior ) );
      const StateT stateF = set2State( species, rotateToFace( face, qfarfield ) );

      // acoustic and convective speeds
      const Real ai = sqrt( stateI.speedOfSound2() );
      const Real af = sqrt( stateF.speedOfSound2() );

      const Real uni = stateI.velocity(0);
      const Real unf = stateF.velocity(0);

   // u+/-a travelling acoustic riemann invariants
      const Real rim = uni - 2.*ai*gam1;
      const Real rip = unf + 2.*af*gam1;

      // boundary values of normal velocity and speed of sound
      const Real ub = 0.5*( rim + rip );
      const Real ab = 0.25*(gam-1.)*( rip - rim );

   // linear riemann invariants (entropy and vorticity) are upwind values

   // upwind state for linear fields
      const StateT stateU = set2State( species, rotateToFace( face, quw ) );

   // entropy
      const Real sb =   stateU.speedOfSound2()
                     /( gam*pow( stateU.density(), gam-1. ) );

      // vorticity
      const auto vb = [&]() -> std::array<Real,nDim-1>
     {
         std::array<Real,nDim-1> v;
         for( int i=0; i<nDim-1; ++i ){ v[i]=stateU.velocity(i+1); }
         return v;
     }();

      // density, pressure at boundary
      const Real rb  = pow( ab*ab/(gam*sb), gam1 );
      const Real pb  = rb*ab*ab/gam;

      PrimVarT qb;

      qb[0] = ub;
      for( int i=0; i<nDim-1; ++i ){ qb[i+1]=vb[i]; }
      qb[nDim] = rb;
      qb[nDim+1] = pb;

//    return qfarfield;
      return rotateFromFace( face, set2Set<SolVarT>( species, qb ) );
  }
