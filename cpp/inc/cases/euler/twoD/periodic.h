
# include <conservationLaws/euler/euler.h>

# include <parallalg/array.h>

# include <utils/concepts.h>

# include <cmath>

/*
 * Initialise to constant velocity with top hat scalar distribution
 *    velocity has angle theta from x axis and magnitude speed
 */
   template<EulerVarSet VarT, floating_point Real>
   par::Array<VarT,2> initialise_scalar_tophat( const Mesh<2,Real>& mesh,
                                                const Species<LawType::Euler,Real>& species,
                                                const Real theta,
                                                const Real speed,
                                                const Real lox,
                                                const Real hix,
                                                const Real loy,
                                                const Real hiy )
  {
      constexpr LawType Law = LawType::Euler;
      constexpr BasisType<Law> PrimBasis = BasisType<Law>::Primitive;

      using PrimVarT = VariableSet<Law,2,PrimBasis,Real>;

      par::Array<VarT,2> q(mesh.cells.shape());

   // background velocity
      const Real u = speed*cos( theta );
      const Real v = speed*sin( theta );

   // background density
      const Real r0 = 1.;

   // background speed of sound
      const Real a = speed/species.minf;

   // background pressure
      const Real p = r0*a*a/species.gamma;

   // cell bounds
      const size_t nx = mesh.cells.shape(0);
      const size_t ny = mesh.cells.shape(1);

      for( size_t i=0; i<nx; i++ )
     {
         for( size_t j=0; j<ny; j++ )
        {
            const par::Idx<2> ij{i,j};

            const Real x = mesh.cells[ij].centre[0];
            const Real y = mesh.cells[ij].centre[1];

         // density
            const Real rx = ((x>lox) and (x<hix)) ? 1 : 0;
            const Real ry = ((y>loy) and (y<hiy)) ? 1 : 0;

            const Real r = r0 + 0.05*rx*ry;

            q[ij] = set2Set<VarT>( species, PrimVarT{{u,v,r,p}} );
        }
     }

      return q;
  }
   
