
# include <ode.h>
# include <controls.h>

# include <spatial/muscl.h>

# include <spatial/boundary/boundaryCondition.h>

# include <timestepping/rungeKutta.h>

# include <limiters/limiter.h>
# include <conservationLaws/scalarAdvection/scalarAdvection.h>
# include <conservationLaws/scalarAdvection/boundaryConditions.h>

# include <solutionField/solutionField.h>

# include <mesh/generate/oneD.h>
# include <mesh/mesh.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <vector>
# include <iostream>
# include <fstream>

# include <cmath>

/*
 * One dimensional periodic test cases for the scalar advection equation
 */

constexpr LawType Law = LawType::ScalarAdvection;
constexpr int nDim = 1;
using BasisT = BasisType<Law>;
using Real = double;

// ------- User Inputs -------

// monochromatic soundwave
constexpr Real amplitude = 1.0;
constexpr unsigned int nwavelengths=1;
constexpr Real velocity = 1.;

// discretisation
constexpr int  nx  = 8;
constexpr int  nt  = 5;
constexpr Real cfl = 0.9;

constexpr BasisT SolutionBasis = BasisT::Conserved;

using Flux = RusanovFlux<Law>;

//using Limiter = Limiters::Cada3;
//using Limiter = Limiters::MonotonizedCentral2;
using Limiter = Limiters::NoLimit1;


// ------ typedefs -------------------

using SolVarSet   = VariableSet<  Law,nDim,SolutionBasis,Real>;
using SolVarDelta = VariableDelta<Law,nDim,SolutionBasis,Real>;
using SolVarGrad  = std::array<SolVarDelta,nDim>;
using SolField    = SolutionField<SolVarSet,nDim>;

using FluxRes = FluxResult<Law,nDim,Real>;

using MeshT = Mesh<nDim,Real>;
using Face  = MeshT::Face;


// ------- i/o -------

   void writeState( std::ofstream& os,
                    const Species<Law,  Real>& species,
                    const State<  Law,1,Real>&   state )
  {
      os << state.velocity(0) << " "
         << state.scalar()    << std::endl;
  }


// ------- run script -------

   int main()
  {
      const par::Shape<nDim> cellShape{nx};

      const ODE::Explicit::RungeKutta<Real> rk = ODE::Explicit::ssp34<Real>();
      const UnsteadyTimeControls<Real> timeControls{.nTimesteps=nt, .cfl=cfl};

   // setup
      const Species<Law,Real> species{};

   // initialise mesh
      const MeshT mesh = make_linspace_mesh<Real>( cellShape, 0,nx );

   // create solution field
      SolField q( cellShape );

   // initialise interior to sinusoid
      const auto initial_conditions = [&]( const MeshT::Cell& c ) -> SolVarSet
     {
         const Real x = c.centre[0]/nx;
         const Real l = nwavelengths*x;
         return { velocity, amplitude*cos(2.*M_PI*l) };
     };
      par::transform( initial_conditions,
                      q.interior,
                      mesh.cells );

   // initialise boundaries
      for( SolField::VarField& v : q.boundary ){ par::fill( v, SolVarSet{velocity,1.} ); }
//    for( SolField::BCType&  bc : q.bcTypes ){ bc = SolField::BCType::Periodic; }
      for( SolField::BCType&  bc : q.bcTypes ){ bc = SolField::BCType::Riemann; }

   // boundary condition set
      const std::tuple boundaryConditions{make_periodic_BCond<Law>(),
                                          make_ghostCell_BCond<Law,BoundaryType<Law>::Riemann>()};

   // high order reconstruction and flux functions
      const auto hoflux = make_muscl_flux<Law>( Limiter{}, Flux{} );

   // integrate forward in time
      integrate( timeControls, rk,
                 hoflux, boundaryConditions,
                 species,
                 mesh, q );

   // write solution to file
     {
      std::ofstream solutionFile( "data/scalarPeriodic1D.dat" );

      if( solutionFile.is_open() )
     {
         solutionFile << std::scientific;
         solutionFile.precision(8);

         par::for_each( // write state to file
                        [&]( const SolVarSet& q0 ) -> void
                           { writeState( solutionFile, species, set2State( species, q0 ) ); },
                        // solution array
                        q.interior
                      );
     }
      else
     {
         std::cout << "cannot open \"data/scalarPeriodic1D.dat\" for writing\n" << std::endl;
         return 1;
     }
      solutionFile.close();
     }

      return 0;
  }

