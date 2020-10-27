
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
      const Real velocity = std::sqrt( sref.velocity2() );
//    const Real velocity = 1.;
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
//    d is height of channel
//    l is length of channel before/after bump
//    w is width of bump
//    h is height of bump
constexpr size_t nl = 48;
constexpr size_t nw = 48;
constexpr size_t nd = 96;

constexpr Real d = 2.0;
constexpr Real l = 3.0;
constexpr Real w = 2.0;
constexpr Real h = 0.100;

constexpr size_t nx = 2*nl+nw;
constexpr size_t ny = nd;
//constexpr size_t nx = 30;
//constexpr size_t ny = 30;

// time discretisation
constexpr size_t nt = 6500;
constexpr Real  cfl = 0.4;

# ifdef _OPENMP
constexpr int nthreads=8;
# endif

// flow conditions
constexpr Real mach     = 0.4;
constexpr Real pressure = 1.;
constexpr Real density  = 1.;
constexpr Real angle = 0.00*M_PI;

constexpr Real speed = mach*std::sqrt( 1.4*pressure/density );
constexpr Real velx  = speed*cos( angle );
constexpr Real vely  = speed*sin( angle );

constexpr BasisT SolutionBasis = BasisT::Primitive;

//using Flux = RoeFlux<Law>;
using Flux = RusanovFlux<Law>;

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


// ------ setup -------------------

   // openMP
# ifdef _OPENMP
      omp_set_num_threads(nthreads);
      std::cout << "OpenMP enabled with: " << nthreads << " threads" << std::endl;
# endif

      const par::DualShape<nDim> cellShape{nx,ny};

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
      const MeshT mesh = make_channelbump_mesh<Real>( {d,l,w,h}, {nd,nl,nw} );

   // initialise solution
      const SolVarSet qref = set2Set<SolVarSet>( species,
                                                 PrimVarSet{velx,
                                                            vely,
                                                            density,
                                                            pressure} );
      SolField q(cellShape);
      par::fill( q.interior, qref );

   // boundary conditions
      const std::tuple boundaryConditions{make_ghostCell_BCond<Law,BCType::Riemann>(),
                                          make_ghostCell_BCond<Law,BCType::Entropy>(),
                                          make_flux_BCond<Law,BCType::InviscidWall>(),
                                          make_periodic_BCond<Law>(),
                                          make_fixed_BCond<Law>()};

      for( SolField::VarField& qb : q.boundary ){ par::fill( qb, qref ); }
      q.bcTypes[0] = BCType::Fixed;
      q.bcTypes[1] = BCType::Fixed;
      q.bcTypes[2] = BCType::InviscidWall;
      q.bcTypes[3] = BCType::InviscidWall;

   // high order reconstruction and flux functions
      const auto hoflux = make_muscl_flux<Law>( Limiter{}, Flux{} );

   // integrate forward in time
      integrate( timeControls, rk,
                 hoflux, boundaryConditions,
                 species,
                 mesh, q );

   // write solution to file
      if( false )
     {
         std::ofstream solutionFile( "data/euler/channelBump2D/result.dat" );
   
         if( solutionFile.is_open() )
        {
            auto writePoint = [&]( const auto& p ) -> void
           {
               solutionFile << p[0] << " "
                            << p[1] << " ";
           };

            solutionFile << std::scientific;
            solutionFile.precision(20);

            const auto sref = set2State( species, qref );

            auto writer = [&] ( const par::DualIdx2 idx,
                                const MeshT::Cell&   c0,
                                const SolVarSet&     qc ) -> void
           {
               if( idx[1]!=0 ){ return; }
               writePoint(c0.centre);
               writeState( solutionFile, species, set2State( species, qc ), sref );
//             if( idx[1]==ny-1 ){ solutionFile << "\n"; }
               return;
           };
   
//          for( size_t i=0; i<q.interior.shape(1); ++i )
            for( size_t i=0; i<1; ++i )
           {
               const auto dcentre =  mesh.cells[{1,i}].centre
                                   - mesh.cells[{0,i}].centre;

               const auto ghost_centre = mesh.cells[{0,i}].centre - 2.*dcentre;
               writePoint(ghost_centre);

               writeState( solutionFile,
                           species,
                           set2State( species,
                                      q.boundary[0][{i,1}] ),
                           sref );
           }
//          for( size_t i=0; i<q.interior.shape(1); ++i )
            for( size_t i=0; i<1; ++i )
           {
               const auto dcentre =  mesh.cells[{1,i}].centre
                                   - mesh.cells[{0,i}].centre;

               const auto ghost_centre = mesh.cells[{0,i}].centre - dcentre;
               writePoint(ghost_centre);

               writeState( solutionFile,
                           species,
                           set2State( species,
                                      q.boundary[0][{i,0}] ),
                           sref );
           }
            par::for_each_idx( writer,
                               mesh.cells,
                               q.interior );
        }
         else
        {
            std::cout << "cannot open \"data/euler/channelBump2D/result.dat\" for writing\n" << std::endl;
            return 1;
        }
         solutionFile.close();
     }

   // write solution to file
      if( true )
     {
         std::ofstream solutionFile( "data/euler/channelBump2D/result.dat" );
   
         if( solutionFile.is_open() )
        {
            auto writePoint = [&]( const auto& p ) -> void
           {
               solutionFile << p[0] << " "
                            << p[1] << " ";
           };

            solutionFile << std::scientific;
            solutionFile.precision(20);

            const auto sref = set2State( species, qref );

            auto writer = [&] ( const par::DualIdx2 idx,
                                const MeshT::Cell&   c0,
                                const SolVarSet&     qc ) -> void
           {
               writePoint(c0.centre);
               writeState( solutionFile, species, set2State( species, qc ), sref );
               if( idx[1]==ny-1 ){ solutionFile << "\n"; }
               return;
           };
   
            if( false ) // write left boundary
           {
            for( size_t i=0; i<q.interior.shape(1); ++i )
           {
               const auto dcentre =  mesh.cells[{1,i}].centre
                                   - mesh.cells[{0,i}].centre;

               const auto ghost_centre = mesh.cells[{0,i}].centre - 2.*dcentre;
               writePoint(ghost_centre);

               writeState( solutionFile,
                           species,
                           set2State( species,
                                      q.boundary[0][{i,1}] ),
                           sref );
               if( i==ny-1 ){ solutionFile << "\n"; }
           }
            for( size_t i=0; i<q.interior.shape(1); ++i )
           {
               const auto dcentre =  mesh.cells[{1,i}].centre
                                   - mesh.cells[{0,i}].centre;

               const auto ghost_centre = mesh.cells[{0,i}].centre - dcentre;
               writePoint(ghost_centre);

               writeState( solutionFile,
                           species,
                           set2State( species,
                                      q.boundary[0][{i,0}] ),
                           sref );
               if( i==ny-1 ){ solutionFile << "\n"; }
           }
           }

            par::for_each_idx( writer,
                               mesh.cells,
                               q.interior );
        }
         else
        {
            std::cout << "cannot open \"data/euler/channelBump2D/result.dat\" for writing\n" << std::endl;
            return 1;
        }
         solutionFile.close();
     }

   // write mesh to file
      if( false )
     {
         std::ofstream meshFile( "data/euler/channelBump2D/mesh.dat" );
   
         if( meshFile.is_open() )
        {
            auto writePoint = [&]( const auto& p ) -> void
           {
               meshFile << p[0] << " "
                        << p[1] << " ";
           };

            meshFile << std::scientific;
            meshFile.precision(8);

            auto writer = [&] ( const par::PrimalIdx2 idx,
                                const MeshT::Node&      p ) -> void
           {
               writePoint(p);
               meshFile << "\n";
               if( idx[1]==ny ){ meshFile << "\n"; }
               return;
           };
   
            par::for_each_idx( writer,
                               mesh.nodes );
        }
         else
        {
            std::cout << "cannot open \"data/euler/channelBump2D/mesh.dat\" for writing\n" << std::endl;
            return 1;
        }
         meshFile.close();
     }

      return 0;
  }

