
# include <ode.h>
# include <controls.h>

# include <spatial/muscl.h>
# include <limiters/limiter.h>

# include <timestepping/rungeKutta.h>

# include <conservationLaws/euler/euler.h>
# include <conservationLaws/euler/boundaryConditions.h>

# include <spatial/boundary/boundaryCondition.h>

# include <solutionField/solutionField.h>

# include <mesh/generate/oneD.h>
# include <mesh/mesh.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <vector>
# include <iostream>
# include <fstream>

# include <cases/euler/oneD/shocktubes.h>

/*
 *    Shock tube problem at very low Mach number
 *       Initial conditions from:
 *          "Improved Flux Formulations for Unsteady Low Mach Number Flows"
 *          Sachdev, Hosangadi, Sankaran, 2012
 *       This case produces two very weak shockwaves, and a very low mach (0.0001) contact discontinuity.
 *
 *       Left state: P = 100028.04 Pa  |   Right state: P = 100000.00 Pa
 *                   U = 0             |                U = 0
 *                   T = 300 K         |                T = 300 K
 *
 *       This test case is a useful case for very low mach number acoustics
 *          Many low mach fluxes intended for aeroacoustics have too small a dissipation on the velocity field (in symmetrised form).
 *          These fluxes will show significant velocity oscillations and small pressure oscillations around the contact discontinuity.
 *          Fluxes with the correct dissipation scaling will produce monotone solutions.
 */

constexpr LawType Law = LawType::Euler;
constexpr int nDim = 1;
using BasisT = BasisType<Law>;

// ------- User Inputs -------

constexpr ShockTube1D problem = ShockTube1D::Sods;

using Real = double;

constexpr int  nx  = 128;
constexpr int  nt  = 480;
constexpr Real cfl = 0.80;

constexpr BasisT SolutionBasis = BasisT::Primitive;

//using Flux = RusanovFlux<Law>;
using Flux = RoeFlux<Law>;

//using Flux = Slau<LowMachScaling::Acoustic,
//                  LowMachScaling::Acoustic>;

using Limiter = Limiters::MonotonizedCentral2;
//using Limiter = Limiters::NoLimit1;


// ------ typedefs -------------------

using SolVarSet = VariableSet<Law,nDim,SolutionBasis,Real>;
using SolField  = SolutionField<SolVarSet,nDim>;
using MeshT     = Mesh<nDim,Real>;


// ------- i/o -------

   void writeState( std::ofstream& os,
                    const Species<Law,  Real>& species,
                    const State<  Law,1,Real>&   state )
  {
      const Real gamma = species.gamma;
      const Real a = sqrt( state.speedOfSound2() );
      const Real r = state.density();
      const Real u = state.velocity(0);
      const Real p = state.pressure();
      const Real s = a*a / ( gamma*pow(r,gamma-1) );

      os << r   << " "
         << u/a << " "
         << u   << " "
         << p   << " "
         << s   << "\n";
  }


// ------- run script -------

   int main()
  {
      const par::DualShape<nDim> cellShape{nx};

      const ODE::Explicit::RungeKutta<Real> rk = ODE::Explicit::ssp34<Real>();
      const UnsteadyTimeControls<Real> timeControls{.nTimesteps=nt, .cfl=cfl};

   // setup
      const Species<Law,Real> species = get_air_species<Real>();

   // initialise mesh
      const MeshT mesh = make_linspace_mesh<Real>( cellShape, 0,nx );

   // initialise solution to left/right states
      SolField q = shocktube_initial_solution<SolVarSet>( problem, species, cellShape );

   // boundary condition set
      const std::tuple boundaryConditions{make_ghostCell_BCond<Law,BoundaryType<Law>::Riemann>()};

   // high order reconstruction and flux functions
      const auto hoflux = make_muscl_flux<Law>( Limiter{}, Flux{} );

   // integrate forward in time
      integrate( timeControls, rk,
                 hoflux, boundaryConditions,
                 species,
                 mesh, q );

   // write solution to file
      std::ofstream solutionFile = []()
     {
         if( problem == ShockTube1D::Sods )
        {
            return std::ofstream( "data/euler/shocktube/sods.dat" );
        }
         else if( problem == ShockTube1D::LowMach )
        {
            return std::ofstream( "data/euler/shocktube/lowmach.dat" );
        }
     }();

      if( solutionFile.is_open() )
     {
         solutionFile << std::scientific;
         solutionFile.precision(8);

         par::for_each( // write state to file
                        [&]( const SolVarSet& q0 ) -> void
                           { writeState( solutionFile, species, set2State( species, q0 ) ); },
                        // solution array
                        q.interior );
     }
      else
     {
         std::cout << "cannot open \"data/euler/shocktube/<case-name>.dat\" for writing\n" << std::endl;
         return 1;
     }
      solutionFile.close();

      return 0;
  }

