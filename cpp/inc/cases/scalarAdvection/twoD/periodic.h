
# include <conservationLaws/scalarAdvection/scalarAdvection.h>

# include <parallalg/array.h>

# include <utils/concepts.h>

# include <cmath>

/*
 * Initialise to constant velocity with top hat scalar distribution
 *    velocity has angle theta from x axis and magnitude speed
 */
   template<ScalarVarSet VarT, floating_point Real>
   par::Array<VarT,2> initialise_scalar_tophat( const Mesh<2,Real>& mesh,
                                                const Species<LawType::ScalarAdvection,Real>& species,
                                                const Real theta,
                                                const Real speed,
                                                const Real lox,
                                                const Real hix,
                                                const Real loy,
                                                const Real hiy )
  {
      constexpr LawType Law = LawType::ScalarAdvection;
      constexpr BasisType<Law> ConsBasis = BasisType<Law>::Conserved;

      using ConsVarT = VariableSet<Law,2,ConsBasis,Real>;

      par::Array<VarT,2> q(mesh.cells.shape());

   // background velocity
      const Real u = speed*cos( theta );
      const Real v = speed*sin( theta );

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

            const Real sx = ((x>lox) and (x<hix)) ? 1 : 0;
            const Real sy = ((y>loy) and (y<hiy)) ? 1 : 0;

            const Real s = sx*sy;

            q[ij] = set2Set<VarT>( species, ConsVarT{{u,v,s}} );
        }
     }

      return q;
  }
   
