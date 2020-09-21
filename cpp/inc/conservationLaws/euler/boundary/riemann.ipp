
# pragma once

# include <conservationLaws/euler/euler.h>

# include <geometry/geometry.h>

# include <utils/concepts.h>

# include <array>
# include <type_traits>

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
                           const SolVarT&                    qboundary )
  {
      using StateT = State<LawType::Euler,nDim,Real>;
      using PrimVarT = VariableSet<LawType::Euler,nDim,EulerBases::Primitive,Real>;

      const Real gam =species.gamma;
      const Real gam1=species.gamma1;

   // direction from interior cell to boundary face
      const geom::Direction<nDim,Real> toBoundary = face.centre - cellinterior.centre;

   // face normal points into or out of domain?
      const bool face_into_domain = dot( toBoundary, face.metric[0] ) < 0;

   // face normal points from left to right
      const SolVarT& qr = face_into_domain ? qinterior : qboundary;
      const SolVarT& ql = face_into_domain ? qboundary : qinterior;

   // left/right/average states aligned with boundary face
      const StateT stateL = set2State(  species, rotateToFace( face, ql ) );
      const StateT stateR = set2State(  species, rotateToFace( face, qr ) );
      const StateT stateA = roeAverage( species, stateL, stateR );

   // boundary mach number
      const Real ma = stateA.velocity(0) / sqrt( stateA.speedOfSound2() );

   // supersonic cases, return fully upwinded value
      if( ma<-1. ){ return qr; }
      if( ma> 1. ){ return ql; }

   // subsonic cases
      const bool flow_to_right = ma > 0

      // acoustic and convective speeds
      const Real al = sqrt( stateL.speedOfSound2() );
      const Real ar = sqrt( stateR.speedOfSound2() );

      const Real unl = stateL.velocity(0);
      const Real unr = stateR.velocity(0);

      // left/right travelling acoustic riemann invariants
      const Real ril = flow_to_right ? (unl + 2.*al*gam1) : (unr + 2.*ar*gam1);
      const Real rir = flow_to_right ? (unr - 2.*ar*gam1) : (unl - 2.*al*gam1);

      // boundary values
         // normal velocity and speed of sound at boundary
      const Real ub = 0.5*( ril + rir );
      const Real ab = 0.25*(gam-1.)*( ril - rir );

      // entropy at boundary is upwind value
      const Real sb = flow_to_right ? al2/( gam*pow( stateL.density(), gam-1. ) )
                                    : ar2/( gam*pow( stateR.density(), gam-1. ) );

      // vorticity at boundary is upwind value
      const auto get_vorticity = []( const StateT& s ) -> std:array<Real,nDim-1>
     {
         std::array<Real,nDim-1> v;
         for( int i=0; i<nDim-1; ++i )
        {
            v[i]=s.velocity(i+1);
        }
         return v;
     };
      const auto vb = flow_to_right ? get_vorticity( stateL )
                                    : get_vorticity( stateR );

      // density, pressure at boundary
      const Real rb  = pow( ab*ab/(gam*sb), gam1 );
      const Real pb  = rb*ab*ab/gam;

      PrimVarT qb;

      qb[0]  = ub;
      for( int i=0; i<nDim-1; ++i ){ qb[i+1]=vb[i]; }
      qb[nDim] = rb;
      qb[nDim+1] = pb;

      return rotateFromFace( face, set2Set<SolVarT>( species, qb ) );
  }

