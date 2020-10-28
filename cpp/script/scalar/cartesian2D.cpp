
# include <ode.h>
# include <controls.h>

# include <spatial/muscl.h>

# include <spatial/boundary/boundaryCondition.h>

# include <timestepping/rungeKutta.h>

# include <limiters/limiter.h>
# include <conservationLaws/scalarAdvection/scalarAdvection.h>
# include <conservationLaws/scalarAdvection/boundaryConditions.h>

# include <solutionField/solutionField.h>

# include <mesh/generate/twoD.h>
# include <mesh/mesh.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <vector>
# include <iostream>
# include <fstream>

# include <omp.h>

# include <cases/scalarAdvection/twoD/periodic.h>

/*
 * Two dimensional periodic test cases
 */

constexpr LawType Law = LawType::ScalarAdvection;
constexpr int nDim = 2;
using BasisT = BasisType<Law>;
using Real = double;

// ------- User Inputs -------

// Convecting velocity - angle and magnitude
constexpr Real theta = 1.20*M_PI;
constexpr Real speed = 1.0;

// lo/hi x/y coordinates of scalar top-hat
constexpr Real th_lx = 0.375;
constexpr Real th_hx = 0.625;

constexpr Real th_ly = 0.375;
constexpr Real th_hy = 0.625;


// discretisation
constexpr int  nx  = 160;
constexpr int  ny  = 80;
constexpr int  nt  = 20;
constexpr Real cfl = 1.5;

constexpr BasisT SolBasis = BasisT::Conserved;

using Flux = RusanovFlux<Law>;

using Limiter = Limiters::NoLimit3;


// ------ typedefs -------------------

using SolVarSet = VariableSet<  Law,nDim,SolBasis,Real>;
using SolField  = SolutionField<SolVarSet,nDim>;
using MeshT     = Mesh<nDim,Real>;


// ------- i/o -------

   void writeState( std::ostream& os,
                    const MeshT::Cell&                                    cell,
                    const Species<LawType::ScalarAdvection,Real>&      species,
                    const State<  LawType::ScalarAdvection,nDim,Real>&   state )
  {
      const Real x = cell.centre[0];
      const Real y = cell.centre[1];

      const Real u = state.velocity(0);
      const Real v = state.velocity(1);
      const Real q = state.scalar();

      os << x << " "
         << y << "  "
         << u << " "
         << v << " "
         << q << std::endl;
  }

// ------- run script -------

   int main()
  {
   // openMP
# ifdef _OPENMP
      omp_set_num_threads(2);
      std::cout << "OpenMP enabled" << std::endl;
# endif

   // setup
      const par::DualShape<nDim> cellShape{nx,nx};

      const ODE::Explicit::RungeKutta<Real> rk = ODE::Explicit::ssp34<Real>();
      const UnsteadyTimeControls<Real> timeControls{.nTimesteps=nt, .cfl=cfl};

      const Species<Law,Real> species{};

   // initialise mesh
      const MeshT mesh = make_linspace_mesh<Real>( cellShape, 0,nx, 0,ny );

   // initialise solution
      SolField q( cellShape );
      q.interior = initialise_scalar_tophat<SolVarSet>( mesh,
                                                        species,
                                                        theta,speed,
                                                        th_lx*nx,th_hx*nx,
                                                        th_ly*ny,th_hy*ny );

//    const SolVarSet qref{cos(theta)*speed,sin(theta)*speed,1.};
//    par::fill( q.interior, qref );

   // initialise boundaries
      for( SolField::VarField& v : q.boundary ){ par::fill( v, q.interior({0,0}) ); }
      q.bcTypes[0] = SolField::BCType::Periodic;
      q.bcTypes[1] = SolField::BCType::Periodic;
      q.bcTypes[2] = SolField::BCType::Riemann;
      q.bcTypes[3] = SolField::BCType::Riemann;

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
      std::ofstream solutionFile("data/scalar/cartesian2D/result.dat");

      if( solutionFile.is_open() )
     {
         solutionFile << std::scientific;
         solutionFile.precision(8);

         for( size_t i=0; i<q.interior.shape(0); ++i )
        {
            for( size_t j=0; j<q.interior.shape(1); ++j )
           {
               writeState( solutionFile, mesh.cells({i,j}), species, set2State( species, q.interior({i,j}) ) );
           }
            solutionFile << std::endl;
        }
     }
      else
     {
         std::cout << "cannot open \"data/scalar/cartesian2D/result.dat\" for writing\n" << std::endl;
         return 1;
     }
      solutionFile.close();
     }

      return 0;
  }

