
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

# include <utils/maths/misc.h>

# include <vector>
# include <iostream>
# include <fstream>

# include <omp.h>

# include <cases/euler/oneD/periodic.h>

/*
 * One dimensional periodic test cases for the Euler equations
 */

constexpr LawType Law = LawType::Euler;
constexpr int nDim = 1;
using BasisT = BasisType<Law>;
using BCType = BoundaryType<Law>;
using Real = double;


// ------- User Inputs -------

constexpr Periodic1D problem = Periodic1D::Other;

// // isolated soundwave
// constexpr Real amplitude = 0.1;
// constexpr unsigned int nwavelengths=1;
// 
// // acoustic - entropy interaction
// constexpr Real centre_a = 0.2;
// constexpr Real centre_s = 0.6;
// 
// constexpr Real width_a = 0.02;
// constexpr Real width_s = 0.04;
// 
// constexpr Real amplitude_a = 0.05;
// constexpr Real amplitude_s = 0.00*mach;

// discretisation
constexpr int  nx  = 24;
constexpr int  nt  = 6000;
constexpr Real cfl = 0.8;

// left/right contact wave
constexpr Real mach      = 0.4;
constexpr Real pressure  = 1.;
constexpr Real density   = 1.;
constexpr Real velocity  = mach*std::sqrt( 1.4*pressure/density );
constexpr Real amplitude = 0.1;

constexpr BasisT SolutionBasis = BasisT::Primitive;

//using Flux = CentralFlux<Law>;
using Flux = RoeFlux<Law>;

//using Flux = Slau<LowMachScaling::Convective,
//                  LowMachScaling::Acoustic>;

//using Limiter = Limiters::Cada3;
//using Limiter = Limiters::MonotonizedCentral2;
using Limiter = Limiters::NoLimit1;


// ------ typedefs -------------------

constexpr BasisT PrimBasis = BasisT::Primitive;

using SolVarSet  = VariableSet<Law,nDim,SolutionBasis,Real>;
using PrimVarSet = VariableSet<Law,nDim,    PrimBasis,Real>;
using SolField   = SolutionField<SolVarSet,nDim>;
using MeshT      = Mesh<nDim,Real>;


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
      const Real sonic = std::sqrt( gamma*pressure/density );
      const Real entropy = sonic*sonic / ( gamma*pow(density,gamma-1) );

      os << r/density  << " "
         << u/velocity << " "
         << p/pressure << " "
         << s/entropy << " "
         << a/sonic << "\n";
  }


// ------- run script -------

   int main()
  {

# ifdef _OPENMP
      omp_set_num_threads(1);
# endif
      const par::DualShape<nDim> cellShape{nx};

      const ODE::Explicit::RungeKutta<Real> rk = ODE::Explicit::ssp11<Real>();
      const UnsteadyTimeControls<Real> timeControls{.nTimesteps=nt, .cfl=cfl};

   // setup
      const Species<Law,Real> species = []() -> Species<Law,Real>
     {
         Species<Law,Real> gas = get_air_species<Real>();
         gas.minf=mach;
         return gas;
     }();

   // initialise mesh
      const MeshT mesh = make_linspace_mesh<Real>( cellShape, 0,nx );

   // initialise solution
      const SolVarSet ql{velocity, (1.+amplitude)*density, pressure};
      const SolVarSet qr{velocity,                density, pressure};

      SolField q(cellShape);
//    for( size_t i=0;    i<nx/2; ++i ){ q.interior({i})=ql; }
//    for( size_t i=nx/2; i<nx;   ++i ){ q.interior({i})=qr; }

      [[maybe_unused]]
      const auto contact = [&]( const MeshT::Cell& cell ) -> SolVarSet
     {
         Real c0 = utils::gaussian( nx/2.,
                                    nx/10.,
                                    cell.centre[0] );
         if( c0<0.001 ){ c0=0.; }
         const Real p = pressure;
         const Real r = density*(1.+amplitude*c0);
         const Real u = velocity;
         return set2Set<SolVarSet>( species, PrimVarSet{{u,r,p}} );
     };
      [[maybe_unused]]
      const auto acoustic = [&]( const MeshT::Cell& cell ) -> SolVarSet
     {
         const Real a0 = std::sqrt( species.gamma*pressure/density );
         const Real c0 = utils::gaussian( nx/2.,
                                          nx/10.,
                                          cell.centre[0] );
         const Real p = pressure*(1.+amplitude*c0);
         const Real r = density*pow( p/pressure, 1./species.gamma );
         const Real a = std::sqrt( species.gamma*p/r );
         const Real u = velocity - 2.*( a0 - a )/(species.gamma-1.);
         return set2Set<SolVarSet>( species, PrimVarSet{{u,r,p}} );
     };

      [[maybe_unused]]
      const auto constant = [&]( const MeshT::Cell& ) -> SolVarSet { return qr; };

      par::transform( constant,
                      q.interior,
                      mesh.cells );

   // boundary conditions
      const std::tuple boundaryConditions{make_ghostCell_BCond<Law,BCType::Riemann>(),
                                          make_ghostCell_BCond<Law,BCType::Entropy>(),
                                          make_periodic_BCond<Law>()};

      const SolVarSet qref{velocity,density,pressure};
      par::fill( q.interior,    qref );
      par::fill( q.boundary[0], qref );
      par::fill( q.boundary[1], qref );
      q.bcTypes[0] = BCType::Riemann;
      q.bcTypes[1] = BCType::Riemann;
//    q.bcTypes[0] = BCType::Periodic;
//    q.bcTypes[1] = BCType::Periodic;
/*
      q.interior = [&]() -> SolField::VarField
     {
         if constexpr ( problem == Periodic1D::Soundwave )
        {
            return soundwave_initial_solution<SolVarSet>( species,
                                                          mesh.cells,
                                                          Real(nx),
                                                          mach,
                                                          amplitude,
                                                          nwavelengths );
        }
         else if constexpr ( problem == Periodic1D::AcousticEntropy )
        {
            return acoustic_entropy_initial_solution<SolVarSet>( species,
                                                                 mesh.cells,
                                                                 Real(nx),
                                                                 mach,
                                                                 centre_a,
                                                                 centre_s,
                                                                 width_a,
                                                                 width_s,
                                                                 amplitude_a,
                                                                 amplitude_s );
        }
     }();
*/

   // high order reconstruction and flux functions
      const auto hoflux = make_muscl_flux<Law>( Limiter{}, Flux{} );

   // integrate forward in time
      integrate( timeControls, rk,
                 hoflux, boundaryConditions,
                 species,
                 mesh, q );

   // write solution to file
     {
         std::ofstream solutionFile = []() -> std::ofstream
        {
            if( problem == Periodic1D::Soundwave )
           {
               return std::ofstream( "data/euler/oneD/soundwave.dat" );
           }
            else if( problem == Periodic1D::AcousticEntropy )
           {
               return std::ofstream( "data/euler/oneD/acoustic_entropy.dat" );
           }
            else
           {
               return std::ofstream( "data/euler/oneD/result.dat" );
           }
        }();
   
         if( solutionFile.is_open() )
        {
            solutionFile << std::scientific;
            solutionFile.precision(12);
   
            writeState( solutionFile, species, set2State( species, q.boundary[0]({1}) ) );
            writeState( solutionFile, species, set2State( species, q.boundary[0]({0}) ) );
            par::for_each( // write state to file
                           [&]( const SolVarSet& q0 ) -> void
                              { writeState( solutionFile, species, set2State( species, q0 ) ); },
                           // solution array
                           q.interior );
            writeState( solutionFile, species, set2State( species, q.boundary[1]({0}) ) );
            writeState( solutionFile, species, set2State( species, q.boundary[1]({1}) ) );
        }
         else
        {
            std::cout << "cannot open \"data/euler/oneD/<case-name>.dat\" for writing\n" << std::endl;
            return 1;
        }
         solutionFile.close();
     }

      return 0;
  }

