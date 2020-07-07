
# pragma once

# include <conservationLaws/euler/euler.h>

# include <parallalg/array.h>

# include <array>

# include <iostream>

// ---------- free-stream boundary conditions ----------

/*
 * calculate freestream boundary values at either end of mesh
 */

   template<int nDim, LawType Law, BasisType<Law> Basis, floating_point Real>
      requires   ImplementedLawType<Law>
              && (nDim==1)
              && ImplementedVarSet< VariableSet<Law,nDim,Basis,Real>>
   std::array<VariableSet<Law,nDim,Basis,Real>,2> 
           updateBCs( const Species<Law,Real>&                          species,
                      const par::Array<geom::Point<nDim,Real>,1>&            nodes,
                      const par::Array<   VariableSet<Law,nDim,Basis,Real>,1>&   q,
                      const std::array<VariableSet<Law,nDim,Basis,Real>,2>& qb0 )
  {
      using SolVarT = VariableSet<Law,nDim,Basis,Real>;

      std::array<SolVarT,2> qbc;

	// left/right boundary cell centres
      const size_t ic0=0;
      const size_t ic1=q.shape(0)-1;

	// left/right boundary nodes
      const size_t in0=0;
      const size_t in1=nodes.shape(0)-1;

   // left face
      qbc[0] = riemannInvariantBC( species,
                                   surface( nodes[{in0}] ),
                                   qb0[0],
                                   q[{ic0}] );

   // right face
      qbc[1] = riemannInvariantBC( species,
                                   surface( nodes[{in1}] ),
                                   q[{ic1}],
                                   qb0[1] );

      return qbc;
  }

/*
 * calculate boundary value from Riemann invariant conditions
 */
   
   template<EulerVarSet SolVarT, floating_point Real>
      requires SameDim<SolVarT,dim_constant<1>>
   SolVarT riemannInvariantBC( const Species<LawType::Euler,Real>& species,
                               const geom::Surface<1,Real>&           face,
                               const SolVarT&                           ql,
                               const SolVarT&                           qr )
  {
      using StateT = State<LawType::Euler,1,Real>;
      using ConservedVarT = VariableSet<LawType::Euler,1,EulerBases::Conserved,Real>;

      const Real gam =species.gamma;
      const Real gam1=species.gamma1;

      const StateT statel = set2State( species, ql );
      const StateT stater = set2State( species, qr );

      const Real al2= statel.speedOfSound2();
      const Real ar2= stater.speedOfSound2();
      const Real al = sqrt( al2 );
      const Real ar = sqrt( ar2 );
      const Real aa = 0.5*( al+ar );

      const Real unl =  projectedVelocity( face.metric[0], statel );
      const Real unr =  projectedVelocity( face.metric[0], stater );
      const Real una = 0.5*( unl+unr );

      const Real ma = una/aa;

   // supersonic cases, return fully upwinded value

      if( ma<-1. ){ return ql; }
      if( ma> 1. ){ return qr; }

   // subsonic cases

      // left/right travelling riemann invariants
      const Real ril = (ma>0) ? (unl + 2.*al*gam1) : (unr + 2.*ar*gam1);
      const Real rir = (ma>0) ? (unr - 2.*ar*gam1) : (unl - 2.*al*gam1);

      // upstream entropy
      const Real sl =  al2/( gam*pow( statel.density(), gam-1. ) );
      const Real sr =  ar2/( gam*pow( stater.density(), gam-1. ) );

      // boundary values
         // velocity, speed of sound, entropy, density, pressure, total energy
      const Real ub = 0.5*( ril + rir );
      const Real ab = 0.25*(gam-1.)*( ril - rir );
      const Real sb = (ma<0) ? sl : sr;
      const Real rb = pow( ab*ab/(gam*sb), gam1 );
      const Real pb = rb*ab*ab/gam;

		const Real rub = rb*ub;
      const Real reb = pb*gam1 + 0.5*rb*ub*ub;

      return set2Set<SolVarT>( species, ConservedVarT{rub,rb,reb} );
  }








