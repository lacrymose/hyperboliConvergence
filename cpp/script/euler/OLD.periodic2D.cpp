
# include <timestepping/rungeKutta.h>

# include <spatial/gradientCalc2D.h>
# include <spatial/residualCalc2D.h>

# include <limiters/limiter.h>
# include <conservationLaws/scalarAdvection/scalarAdvection.h>
# include <conservationLaws/euler/euler.h>

# include <mesh/generate/twoD.h>
# include <mesh/mesh.h>

# include <controls.h>

# include <ode.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <vector>
# include <iostream>
# include <fstream>

# include <cases/scalarAdvection/twoD/periodic.h>
# include <cases/euler/twoD/periodic.h>

/*
 * Two dimensional periodic test cases
 */

//constexpr LawType Law = LawType::ScalarAdvection;
constexpr LawType Law = LawType::Euler;
constexpr int nDim = 2;
using BasisT = BasisType<Law>;
using Real = double;

// ------- User Inputs -------

// Convecting velocity - angle and magnitude
constexpr Real theta = 0.25*M_PI;
constexpr Real speed = 1.0;

// lo/hi x/y coordinates of scalar top-hat
constexpr Real th_lx = 0.375;
constexpr Real th_hx = 0.625;

constexpr Real th_ly = 0.375;
constexpr Real th_hy = 0.625;


// discretisation
constexpr int  nx  = 32;
constexpr int  nt  = 128;
constexpr Real cfl = 1.0;

constexpr BasisT SolBasis = BasisT::Primitive;

//using Flux = CentralFlux<Law>;
//using Flux = RoeFlux<Law>;
//using Flux = RusanovFlux<Law>;

using Flux = Slau<LowMachScaling::Convective,
                  LowMachScaling::Acoustic>;

using Limiter = Limiters::NoLimit1;


// ------ typedefs -------------------

using SolVarSet  = VariableSet<  Law,nDim,SolBasis,Real>;
using SolVarDel  = VariableDelta<Law,nDim,SolBasis,Real>;
using SolVarGrad = std::array<SolVarDel,nDim>;

using FluxRes = FluxResult<Law,nDim,Real>;

using MeshT = Mesh<nDim,Real>;
using Face  = MeshT::Face;


// ------- i/o -------

   void writeState( std::ofstream& os,
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

   void writeState( std::ofstream& os,
                    const MeshT::Cell&                          cell,
                    const Species<LawType::Euler,Real>&      species,
                    const State<  LawType::Euler,nDim,Real>&   state )
  {
      const Real x = cell.centre[0];
      const Real y = cell.centre[1];

      const Real u = state.velocity(0);
      const Real v = state.velocity(1);
      const Real p = state.pressure();

      const Real gam = species.gamma;
      const Real a2= state.speedOfSound2();
      const Real r = state.density();
      const Real s = a2 / ( gam*pow(r,gam-1.) );

      os << x << " "
         << y << "  "
         << u << " "
         << v << " "
         << s << " "
         << r << " "
         << p << std::endl;
  }


// ------- run script -------

   int main()
  {
//    const ODE::Explicit::RungeKutta<Real>      rk = IO::read_exRungeKutta<Real>(         "data/twoD/euler/periodic2D/rungekutta.dat" );
//    const UnsteadyTimeControls<Real> timeControls = IO::read_UnsteadyTimeControls<Real>( "data/twoD/euler/periodic2D/time.dat" );
//    const Species<Law,Real>               species = IO::read_Species<Law,Real>(          "data/twoD/euler/periodic2D/species.dat" );

      const par::Shape<nDim> cellShape{nx,nx};

      const ODE::Explicit::RungeKutta<Real> rk = ODE::Explicit::ssp33<Real>();
      const UnsteadyTimeControls<Real> timeControls{nt,cfl};

   // setup
      const Species<Law,Real> species=[](){ auto s = get_air_species<Real>(); s.minf=0.1; return s; }();
      const Flux                 flux{};
      const Limiter           limiter{};


   // initialise mesh
      const MeshT mesh = make_linspace_mesh<Real>( cellShape, 0,nx, 0,nx );

   // initialise solution
//    par::Array<SolVarSet,nDim> q = initialise_scalar_tophat<SolVarSet>( mesh,
      par::Array<SolVarSet,nDim> q = initialise_scalar_gaussian<SolVarSet>( mesh,
                                                                          species,
                                                                          theta,speed,
                                                                          th_lx*nx,th_hx*nx,
                                                                          th_ly*nx,th_hy*nx );

   // high order reconstruction and flux functions
      const auto hoflux = [&species, &limiter, &flux]
                          ( const Face&       face,
                            const SolVarDel& gradl,
                            const SolVarDel& gradr,
                            const SolVarSet&    ql,
                            const SolVarSet&    qr ) -> FluxRes
     {
      // central and l/r biased differences in solution basis
         const SolVarDel dqc = qr - ql;
         const SolVarDel dql = gradl - dqc;
         const SolVarDel dqr = gradr - dqc;

      // limited differences in solution basis aligned with background coordinate system
         const SolVarDel slopel = limiter( dqc, dql );
         const SolVarDel sloper = limiter( dqc, dqr );

      // inviscid flux
         return flux( species, face, ql+0.5*slopel,
                                     qr-0.5*sloper );
     };

   // gradient array
      par::Array<SolVarGrad,nDim>   dq(cellShape);

   // spatial residual evaluation
      const auto rhs = [&dq,&hoflux]
                       ( const MeshT&                    msh,
                         const Real                     time,
                         const par::Array<SolVarSet,nDim>& v )
     {
         dq = gradientCalcP( msh.cells, v );
         return residualCalc( msh, hoflux, v, dq );
     };

      integrate( timeControls, rk,
                 species, mesh, q,
                 rhs );

   // write solution to file
     {
      std::ofstream solutionFile = []()
     {
         if( Law==LawType::ScalarAdvection ){ return std::ofstream("data/twoD/periodic/scalarAdvection/result.dat"); }
         else if( Law==LawType::Euler ){ return std::ofstream("data/twoD/periodic/euler/result.dat"); }
     }();

      if( solutionFile.is_open() )
     {
         solutionFile << std::scientific;
         solutionFile.precision(12);

         for( size_t i=0; i<q.shape(0); i++ )
        {
            for( size_t j=0; j<q.shape(1); j++ )
           {
               writeState( solutionFile, mesh.cells({i,j}), species, set2State( species, q({i,j}) ) );
           }
            solutionFile << std::endl;
        }
     }
      else
     {
         std::cout << "cannot open \"data/twoD/periodic/scalarAdvection/result.dat\" for writing\n" << std::endl;
         return 1;
     }
      solutionFile.close();
     }

      return 0;
  }

