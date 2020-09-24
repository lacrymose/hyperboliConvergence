
# include <solutionField/solutionField.h>

# include <conservationLaws/euler/euler.h>

# include <geometry/geometry.h>

# include <parallalg/array.h>

# include <cassert>

/*
 * Initialisation routines for 1 dimensional shock tube problems for the Euler equations
 *
 *
 *    An enum class is used to tag each problem, and all functions are templated by the tag.
 *    This means that same script can be used for all problems, just changing which tag is used.
 *
 *    For each problem, several functions are provided. In the descriptions below:
 *       VarSet is a 1D, Euler claw::VariableSet
 *       Species is an Euler claw::Species
 *       Cell is a 1D geom::Volume
 *
 *    Return the initial left/right state for the problem TAG in primitive variables with floating point type Real
 *
 *       shocktube_initial_left< VarSet>( ShockTube1D, Species& ) -> VarSet;
 *       shocktube_initial_right<VarSet>( ShockTube1D, Species& ) -> VarSet;
 *
 *    Return an full 1D array of initial conditions over a given array of cells.
 *    Cell array must have an even length (otherwise diaphragm will sit inside a cell)
 *
 *       shocktube_initial_solution<VarSet>( ShockTube1D, Species&, par::Array<Cell,1>& ) -> par::Array<VarSet,1>
 *    
 *    !!! NOT IMPLEMENTED YET !!!
 *    Exact solution at time t
 *
 *       shocktube_exact_solution<VarSet>( ShockTube1D, Species&, par::Array<Cell,1>&, Real t ) -> par::Array<VarSet,1>
 *    !!! NOT IMPLEMENTED YET !!!
 */

   enum struct ShockTube1D
  {
      Sods,
      LowMach,
//    Toro1,
//    Toro2,
//    Toro3,
//    Toro4,
//    Toro5,
//    Toro6,
//    Toro7,
  };

   template<EulerVarSet VarT, floating_point Real>
      requires SameFPType<VarT,Real>
   VarT shocktube_initial_left( const ShockTube1D                   problem,
                                const Species<LawType::Euler,Real>& species )
  {
      using PrimVarT = VariableSet<LawType::Euler,1,EulerBases::Primitive,Real>;

      switch( problem )
     {
         case( ShockTube1D::Sods ) :
        {
            constexpr Real ul = 0.;
            constexpr Real rl = 1.;
            constexpr Real pl = 1.;
            return set2Set<VarT>( species, PrimVarT{{ul,rl,pl}} );
        }

         case( ShockTube1D::LowMach ) :
        {
            constexpr Real ul = 0.;
            constexpr Real tl = 300.;
            constexpr Real pl = 100028.04;
            const Real rl = pl / ( species.R * tl );
            return set2Set<VarT>( species, PrimVarT{{ul,rl,pl}} );
        }
     }
      return set2Set<VarT>( species, PrimVarT{{0.,0.,0.}} );
  }


   template<EulerVarSet VarT, floating_point Real>
      requires SameFPType<VarT,Real>
   VarT shocktube_initial_right( const ShockTube1D                   problem,
                                 const Species<LawType::Euler,Real>& species )
  {
      using PrimVarT = VariableSet<LawType::Euler,1,EulerBases::Primitive,Real>;

      switch( problem )
     {
         case( ShockTube1D::Sods ) :
        {
            constexpr Real ur = 0.;
            constexpr Real rr = 0.125;
            constexpr Real pr = 0.1;
            return set2Set<VarT>( species, PrimVarT{{ur,rr,pr}} );
        }

         case( ShockTube1D::LowMach ) :
        {
            constexpr Real ur = 0.;
            constexpr Real tr = 300.;
            constexpr Real pr = 100000.;
            const Real rr = pr / ( species.R * tr );
            return set2Set<VarT>( species, PrimVarT{{ur,rr,pr}} );
        }
     }
      return set2Set<VarT>( species, PrimVarT{{0.,0.,0.}} );
  }

   template<EulerVarSet VarT, floating_point Real>
      requires SameFPType<VarT,Real>
   std::array<VarT,2> shocktube_initial_leftright( const ShockTube1D                   problem,
                                                   const Species<LawType::Euler,Real>& species )
  {
      return { shocktube_initial_left< VarT>( problem, species ),
               shocktube_initial_right<VarT>( problem, species ) };
  }

   template<EulerVarSet VarT, floating_point Real>
      requires SameFPType<VarT,Real>
   SolutionField<VarT,1> shocktube_initial_solution( const ShockTube1D                    problem,
                                                  const Species<LawType::Euler,Real>&  species,
                                                  const par::Shape<1>                meshShape )
  {
      const VarT ql = shocktube_initial_left< VarT>( problem, species );
      const VarT qr = shocktube_initial_right<VarT>( problem, species );

      const size_t nx = meshShape[0];
      assert( nx%2 == 0 );

      SolutionField<VarT,1> q(meshShape);

      for( size_t i=0;    i<nx/2; ++i ){ q.interior[{i}] = ql; }
      for( size_t i=nx/2; i<nx;   ++i ){ q.interior[{i}] = qr; }

      for( EulerBCs& bc : q.bcTypes ){ bc=EulerBCs::Riemann; }
      par::fill( q.boundary[0], ql );
      par::fill( q.boundary[1], qr );

      return q;
  }






