
# include <spatial/residualCalc.h>
# include <spatial/spectralRadius.h>
# include <spatial/eulerForwardUpdate.h>

# include <conservationLaws/scalarAdvection/scalarAdvection.h>
# include <conservationLaws/euler/euler.h>

# include <geometry/geometry.h>

# include <vector>
# include <iostream>
# include <fstream>

# include <cassert>


// ------- Inputs ------- 
using Real = float;

constexpr LawType Law = LawType::Euler;
constexpr int nDim = 1;

constexpr int  nx  = 2048;
constexpr int  nt  = 2048;
constexpr Real cfl = 0.5;

using BasisT = BasisType<Law>;
//constexpr BasisT SolutionBasis = BasisT::Conserved;
constexpr BasisT SolutionBasis = BasisT::Viscous;

constexpr bool print=false;

// -----------------------------------

constexpr BasisT ConservedBasis = BasisT::Conserved;

using SolutionVarSet    = VariableSet<  Law,nDim,SolutionBasis, Real>;
using SolutionVarDelta  = VariableDelta<Law,nDim,SolutionBasis, Real>;
using ConservedVarDelta = VariableDelta<Law,nDim,ConservedBasis,Real>;

using FluxRes = FluxResult<Law,nDim,Real>;

using Flux = RusanovFlux<Law>;

using Point   = geom::Point<  nDim,Real>;
using Surface = geom::Surface<nDim,Real>;
using Volume  = geom::Volume< nDim,Real>;

  
# include "../script/riemannProblem.h"

   int main()
  {

      std::cout << std::scientific;
      std::cout.precision(3);

   // setup
      const Species<Law,Real> species = get_species( ConservedBasis );
      const Flux flux{};

   // solutions arrays
      std::vector<SolutionVarSet>  q(nx);
      std::vector<SolutionVarSet> q1(nx);
      std::vector<FluxRes>       res(nx);

   // initialise mesh
      const std::vector<Point> nodes = []()
     {
         std::vector<Point> nds(nx+1);
         int i=0;
         for( Point& p : nds ){ p = Point{ Real(i++) }; };
         return nds;
     }();
      const std::vector<Volume> cells = geom::dual( nodes );

   // initialise solution to left/right states
      const SolutionVarSet ql0 = set2Set<SolutionVarSet>( species, initialLeft(  SolutionBasis ) );
      const SolutionVarSet qr0 = set2Set<SolutionVarSet>( species, initialRight( SolutionBasis ) );

      std::cout << "left  state: " << ql0 << std::endl;
      std::cout << "right state: " << qr0 << std::endl;

      for( int i=0;    i<nx/2; i++ ){ q[i]=ql0; }
      for( int i=nx/2; i<nx;   i++ ){ q[i]=qr0; }

   // integrate forward in time
      for( int t=0; t<nt; t++ )
     {
      // accumulate flux residual
         residualCalc( cells, nodes, flux, species, ql0, q, res );

      // calculate maximum stable timestep
         const Real lmax = spectralRadius( cells, res );

      // integrate cell residuals forward by dt and average over cell volume
         eulerForwardUpdate( cells, species, cfl, lmax, res, q, q1 );
         std::swap( q,q1 );
     }

      std::ofstream solutionFile( "data/testSA.dat" );
      if( solutionFile.is_open() )
     {
         for( const SolutionVarSet& q0 : q )
        {
            writeState( solutionFile, set2State( species, q0 ) );
        }
     }
      else
     {
         std::cout << "cannot open \"testSA.dat for writing\"\n" << std::endl;
         return 1;
     }
      solutionFile.close();

      return 0;
  }

