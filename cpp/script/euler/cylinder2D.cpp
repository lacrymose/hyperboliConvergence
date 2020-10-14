
# include <ode.h>
# include <controls.h>

# include <spatial/muscl.h>
# include <limiters/limiter.h>

# include <timestepping/rungeKutta.h>

# include <conservationLaws/euler/euler.h>
# include <conservationLaws/euler/boundaryConditions.h>

# include <spatial/boundary/boundaryCondition.h>

# include <solutionField/solutionField.h>

# include <mesh/generate/twoD.h>
# include <mesh/mesh.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <utils/maths/misc.h>

# include <vector>
# include <iostream>
# include <fstream>

# include <omp.h>

# include <cases/euler/oneD/periodic.h>

// ------- i/o -------

   template<floating_point Real>
   void writeState( std::ofstream& os,
                    const Species<LawType::Euler,  Real>& species,
                    const State<  LawType::Euler,2,Real>&   state,
                    const State<  LawType::Euler,2,Real>&    sref )
  {
      const Real gamma = species.gamma;

      const Real a = sqrt( state.speedOfSound2() );
      const Real r = state.density();
      const Real u = state.velocity(0);
      const Real v = state.velocity(1);
      const Real p = state.pressure();
      const Real s = a*a / ( gamma*pow(r,gamma-1) );

      const Real density  = sref.density();
      const Real pressure = sref.pressure();
//    const Real velocity = std::sqrt( sref.velocity2() );
      const Real velocity = 1.;
      const Real sonic    = std::sqrt( sref.speedOfSound2() );
      const Real entropy  = sonic*sonic / ( gamma*pow(density,gamma-1) );

      os << r/density  << " "
         << u/velocity << " "
         << v/velocity << " "
         << p/pressure << " "
         << s/entropy  << " "
         << a/sonic    << "\n";
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

// space discretisation
//    rc is radius of cylinder
//    rb is radius of boundary
//    th is theta direction
//    h is height of bump
constexpr size_t nr  = 48;
constexpr size_t nth = 48;

constexpr Real rc =  1.0;
constexpr Real rb = 10.0;

// time discretisation
constexpr size_t nt = 104;
constexpr Real  cfl = 0.4;

# ifdef _OPENMP
constexpr int nthreads=4;
# endif

// flow conditions
constexpr Real mach     = 0.5;
constexpr Real pressure = 1.;
constexpr Real density  = 1.;

constexpr Real velocity  = mach*std::sqrt( 1.4*pressure/density );

constexpr BasisT SolutionBasis = BasisT::Primitive;

//using Flux = RoeFlux<Law>;
//using Flux = RusanovFlux<Law>;

using Flux = Slau<LowMachScaling::Convective,
                  LowMachScaling::Acoustic>;

//using Limiter = Limiters::VanAlbada2;
//using Limiter = Limiters::MonotonizedCentral2;
using Limiter = Limiters::NoLimit1;


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

      const ODE::Explicit::RungeKutta<Real> rk = ODE::Explicit::ssp11<Real>();
      const UnsteadyTimeControls<Real> timeControls{.nTimesteps=nt, .cfl=cfl};

      const Species<Law,Real> species = []() -> Species<Law,Real>
     {
         Species<Law,Real> gas = get_air_species<Real>();
         gas.minf=mach;
         return gas;
     }();

   // initialise mesh
//    const MeshT mesh = make_linspace_mesh<Real>( cellShape, -(l+0.5*w),(l+0.5*w), 0,d );
      const MeshT mesh = make_cylinder_mesh<Real>( {rc,rb}, {nr,nth} );
      const par::Shape<nDim> cellShape(mesh.cells.shape());

   // initialise solution
//    const SolVarSet q0   = set2Set<SolVarSet>( species, PrimVarSet{Real(0.),Real(0.),density,pressure} );
//    const SolVarSet qref = set2Set<SolVarSet>( species, PrimVarSet{Real(0.),Real(0.),density,pressure} );
      const SolVarSet qref = set2Set<SolVarSet>( species, PrimVarSet{velocity, Real(0.), density, pressure} );
//    const SolVarSet qref = set2Set<SolVarSet>( species, PrimVarSet{Real(0.), velocity, density, pressure} );
//    const SolVarSet qref = set2Set<SolVarSet>( species, PrimVarSet{velocity, velocity, density, pressure} );
      SolField q(cellShape);
      par::fill( q.interior, qref );

   // boundary conditions
      const std::tuple boundaryConditions{make_ghostCell_BCond<Law,BCType::Riemann>(),
                                          make_ghostCell_BCond<Law,BCType::Entropy>(),
                                          make_flux_BCond<Law,BCType::InviscidWall>(),
                                          make_periodic_BCond<Law>()};

      for( SolField::VarField& qb : q.boundary ){ par::fill( qb, qref ); }
//    par::fill( q.boundary[0], qref );
//    par::fill( q.boundary[1], qref );
      q.bcTypes[0] = BCType::InviscidWall;
      q.bcTypes[1] = BCType::Entropy;
      q.bcTypes[2] = BCType::Periodic;
      q.bcTypes[3] = BCType::Periodic;

   // high order reconstruction and flux functions
      const auto hoflux = make_muscl_flux<Law>( Limiter{}, Flux{} );

   // integrate forward in time
      integrate( timeControls, rk,
                 hoflux, boundaryConditions,
                 species,
                 mesh, q );

   // write solution to file
      if( true )
     {
         std::ofstream solutionFile( "data/euler/cylinder2D/result.dat" );
   
         if( solutionFile.is_open() )
        {
            const Real dr  = ( rb - rc ) / nr;
            const Real dth = ( 2.*M_PI ) / nth;

            auto writeCylindricalCell = [&]( const par::Idx<2>& idx ) -> void
           {
               const Real theta = (idx[1]+0.5)*dth;
               const Real z = 0.;
               const Real r = rc + (idx[0]+0.5)*dr;
               solutionFile << theta << " "
                            << z     << " "
                            << r     << " ";
           };

            solutionFile << std::scientific;
            solutionFile.precision(8);

            const auto sref = set2State( species, qref );

            auto writer = [&] ( const par::Idx<2> idx,
                                const MeshT::Cell& c0,
                                const SolVarSet&   qc ) -> void
           {
               writeCylindricalCell(idx);
               writeState( solutionFile, species, set2State( species, qc ), sref );
               if( idx[1]==mesh.cells.shape(1)-1 ){ solutionFile << "\n"; }
               return;
           };
   
            par::for_each_idx( writer,
                               mesh.cells,
                               q.interior );
        }
         else
        {
            std::cout << "cannot open \"data/euler/cylinder2D/result.dat\" for writing\n" << std::endl;
            return 1;
        }
         solutionFile.close();
     }

      if( true )
     {
         std::ofstream meshFile( "data/euler/cylinder2D/mesh.dat" );
   
         if( meshFile.is_open() )
        {
            auto writePoint = [&]( const auto& p ) -> void
           {
               meshFile << p[0] << " "
                        << p[1] << " ";
           };

            meshFile << std::scientific;
            meshFile.precision(8);

            auto writer = [&] ( const par::Idx<2> idx,
                                const MeshT::Node&  p ) -> void
           {
               writePoint(p);
               meshFile << "\n";
               if( idx[1]==mesh.cells.shape(1) ){ meshFile << "\n"; }
               return;
           };
   
            par::for_each_idx( writer,
                               mesh.nodes );
        }
         else
        {
            std::cout << "cannot open \"data/euler/cylinder2D/mesh.dat\" for writing\n" << std::endl;
            return 1;
        }
         meshFile.close();
     }

      return 0;
  }

