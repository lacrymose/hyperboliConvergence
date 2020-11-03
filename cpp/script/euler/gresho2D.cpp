
# include <ode.h>
# include <controls.h>

# include <spatial/muscl.h>
# include <limiters/limiter.h>

# include <timestepping/rungeKutta.h>

# include <spatial/boundary/boundaryCondition.h>

# include <solutionField/solutionField.h>

# include <conservationLaws/euler/euler.h>
# include <conservationLaws/euler/boundaryConditions.h>

# include <mesh/generate/twoD.h>
# include <mesh/mesh.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>
# include <parallalg/parallalg.h>

# include <utils/maths/misc.h>

# include <omp.h>

# include <vector>
# include <iostream>
# include <fstream>

# include <cmath>

// ------- i/o -------

   template<floating_point Real>
   void writeState( std::ofstream& os,
                    const Species<LawType::Euler,  Real>& species,
                    const State<  LawType::Euler,2,Real>&   state,
                    const State<  LawType::Euler,2,Real>&    sref )
  {
      const Real gamma = species.gamma;

      const Real a  = sqrt( state.speedOfSound2() );
      const Real a0 = sqrt(  sref.speedOfSound2() );

      const Real r  = state.density();
      const Real r0 =  sref.density();

      const Real p  = state.pressure();
      const Real p0 =  sref.pressure();

      const Real u  = state.velocity(0);
      const Real v  = state.velocity(1);

      const Real speed  = sqrt( state.velocity2() );
      const Real speed0 = 1.;

      const Real s  =  a*a  / ( gamma*pow(r, gamma-1) );
      const Real s0 = a0*a0 / ( gamma*pow(r0,gamma-1) );

      os << r/r0      << " "
         << u/speed0  << " "
         << v/speed0  << " "
         << (p-p0)    << " "
         << s/s0      << " "
         << a/a0      << " "
         << speed/a;
  }


/*
 * Two dimensional channel with sinusoidal bump test case for the Euler equations
 */

// ------- run script -------

   int main()
  {

constexpr LawType Law = LawType::Euler;
constexpr int nDim = 2;
using BasisT = BasisType<Law>;
using BCType = BoundaryType<Law>;
using Real = double;

// ------- User Inputs -------

# ifndef _OPENMP
constexpr auto policy = par::execution::seq;
# else
constexpr auto policy = par::execution::seq;
constexpr int nthreads=1;
# endif

// space discretisation
constexpr size_t nx = 256;

constexpr Real h = 1.0;

// time discretisation
constexpr size_t nt = 10;
constexpr Real  cfl = 1.6;

// variable flow conditions
constexpr Real mach = 3.e-1;
constexpr Real vel_inf = 0.1;

// fixed flow conditions
constexpr Real density = 1.;

constexpr BasisT SolutionBasis = BasisT::Primitive;

//using Flux = RoeFlux<Law>;
using Flux = RusanovFlux<Law>;
//using Flux = CentralFlux<Law>;

//using Flux = RoeUnprecWS;

//using Flux = Slau<LowMachScaling::Convective,
//                  LowMachScaling::Acoustic>;

//using Flux = AusmPlusUP<LowMachScaling::Convective,
//                        LowMachScaling::Convective>;

//using Limiter = Limiters::Cada3;
//using Limiter = Limiters::MonotonizedCentral2;
using Limiter = Limiters::NoLimit3;


// ------ typedefs -------------------

constexpr BasisT PrimBasis = BasisT::Primitive;

using SolVarSet  = VariableSet<Law,nDim,SolutionBasis,Real>;
using PrimVarSet = VariableSet<Law,nDim,    PrimBasis,Real>;
using SolField   = SolutionField<SolVarSet,nDim>;
using MeshT      = Mesh<nDim,Real>;


// ------ setup -------------------

   // openMP
# ifdef _OPENMP
      omp_set_num_threads(nthreads);
      std::cout << "OpenMP enabled with: " << nthreads << " threads" << std::endl;
# endif

      const par::DualShape<nDim> cellShape{nx,nx};

      const ODE::Explicit::RungeKutta<Real> rk = ODE::Explicit::ssp34<Real>();
      const UnsteadyTimeControls<Real> timeControls{.nTimesteps=nt, .cfl=cfl};

      const Species<Law,Real> species = [&]() -> Species<Law,Real>
     {
         Species<Law,Real> gas = get_air_species<Real>();
         gas.minf=mach;
         return gas;
     }();

   // initialise mesh
      const MeshT mesh = make_linspace_mesh<Real>( cellShape, -h,h, -h,h );

   // initialise solution
      SolField q(cellShape);

      // reference variable set
      const Real pressure0 = density/( species.gamma*mach*mach );
      const Real pressure = pressure0 - 2. + 4.*log(2.);
      const Real zero{0};
      const SolVarSet qref = set2Set<SolVarSet>( species, PrimVarSet{zero,zero,density,pressure} );

//    std::cout << "speed of sound: " << sqrt( set2State( species, qref ).speedOfSound2() ) << "\n";

      const auto gresho_init = [&]( const MeshT::Cell& cell ) -> SolVarSet
     {
         constexpr Real ln2  = std::log( 2.0 );
         constexpr Real ln02 = std::log( 0.2 );

         const Real x = cell.centre[0];
         const Real y = cell.centre[1];
         const Real r = sqrt( x*x + y*y );

         Real u_theta  = 0.;
         Real pressure = pressure0 - 2. + 4.*ln2;

         if( r < 0.2 )
        {
            u_theta = 5.*r;
            pressure = pressure0 + 12.5*r*r;
        }
         else if( r < 0.4 )
        {
            u_theta = 2.-5.*r;
            pressure = pressure0 + 12.5*r*r + 4.*( 1. - 5.*r - ln02 + std::log( r ) );
        }

         const Real sin_theta = y/r;
         const Real cos_theta = x/r;

         const Real ux = vel_inf - u_theta*sin_theta;
         const Real uy =           u_theta*cos_theta;

         return set2Set<SolVarSet>( species, PrimVarSet{ux,uy,density,pressure} );
     };

      par::transform( gresho_init,
                      q.interior,
                      mesh.cells );

   // boundary conditions
      const std::tuple boundaryConditions{make_periodic_BCond<Law>()};

      for( auto& qb : q.boundary ){ par::fill( qb, qref ); }
      for( auto& bc : q.bcTypes  ){ bc = BCType::Periodic; }

   // high order reconstruction and flux functions
      const auto hoflux = make_muscl_flux<Law>( Limiter{}, Flux{} );

   // integrate forward in time
      integrate( policy, timeControls, rk,
                 hoflux, boundaryConditions,
                 species,
                 mesh, q );

   // write solution to file
      if( true )
     {
         std::ofstream solutionFile( "data/euler/gresho2D/result.dat" );
   
         if( solutionFile.is_open() )
        {
            auto writePoint = [&]( const MeshT::Node& p ) -> void
           {
               solutionFile << p[0] << " "
                            << p[1] << " ";
           };

            solutionFile << std::scientific;
            solutionFile.precision(16);

            const auto sref = set2State( species, qref );

            auto writer = [&]( const par::DualIdx2 idx,
                               const MeshT::Cell&   c0,
                               const SolVarSet&     qc ) -> void
           {
               writePoint(c0.centre);
               writeState( solutionFile, species, set2State( species, qc ), sref );
               solutionFile << "\n";
               if( idx[1]==nx-1 ){ solutionFile << "\n"; }
               return;
           };
   
            par::for_each_idx( writer,
                               mesh.cells,
                               q.interior );
        }
         else
        {
            std::cout << "cannot open \"data/euler/gresho2D/result.dat\" for writing\n" << std::endl;
            return 1;
        }
         solutionFile.close();
     }

   // write mesh to file
      if( false )
     {
         std::ofstream meshFile( "data/euler/gresho2D/mesh.dat" );
   
         if( meshFile.is_open() )
        {
            auto writePoint = [&]( const MeshT::Node& p ) -> void
           {
               meshFile << p[0] << " "
                        << p[1];
           };

            meshFile << std::scientific;
            meshFile.precision(8);

            auto writer = [&] ( const par::PrimalIdx2& idx,
                                const MeshT::Node&       p ) -> void
           {
               writePoint(p);
               meshFile << "\n";
               if( idx[1]==nx ){ meshFile << "\n"; }
               return;
           };
   
            par::for_each_idx( writer,
                               mesh.nodes );
        }
         else
        {
            std::cout << "cannot open \"data/euler/gresho2D/mesh.dat\" for writing\n" << std::endl;
            return 1;
        }
         meshFile.close();
     }

      return 0;
  }

