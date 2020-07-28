
# include <ode.h>

# include <spatial/gradientCalc2D.h>
# include <spatial/residualCalc2D.h>
# include <spatial/spectralRadius.h>
# include <spatial/rungeKuttaAccumulation.h>
# include <spatial/eulerForwardUpdate.h>

# include <limiters/limiter.h>
# include <conservationLaws/scalarAdvection/scalarAdvection.h>
# include <conservationLaws/euler/euler.h>

# include <mesh/generate/twoD.h>
# include <mesh/mesh.h>

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
constexpr int  nt  = 640;
constexpr Real cfl = 1.0;

constexpr BasisT SolBasis = BasisT::Conserved;

//using Flux = CentralFlux<Law>;
using Flux = RoeFlux<Law>;
//using Flux = RusanovFlux<Law>;

using Limiter = Limiters::NoLimit3;


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
      const par::Shape<nDim> cellShape{nx,nx};

      const ODE::Explicit::RungeKutta<Real> rk = ODE::Explicit::ssp33<Real>();

      std::cout << std::scientific;
      std::cout.precision(8);

   // setup
//    const Species<Law,Real> species{};
      const Species<Law,Real> species=[](){ auto s = get_air_species<Real>(); s.minf=0.3; return s; }();
      const Flux                 flux{};
      const Limiter           limiter{};

   // initialise mesh
      const MeshT mesh = make_linspace_mesh<Real>( cellShape, 0,nx, 0,nx );

   // initialise solution
      par::Array<SolVarSet,nDim> q = initialise_scalar_tophat<SolVarSet>( mesh,
                                                                          species,
                                                                          theta,speed,
                                                                          th_lx*nx,th_hx*nx,
                                                                          th_ly*nx,th_hy*nx );

      par::Array<SolVarSet,nDim> q1 = par::copy( q );
      par::Array<SolVarSet,nDim> q2 = par::copy( q );

   // gradient array
      par::Array<SolVarGrad,nDim>   dq(cellShape);
      par::fill( dq, SolVarGrad{} );

   // residual evaluated and accumulated at each runge-kutta stage
      // MDArray has no copy constructor so must construct in-place
      par::Array<FluxRes,nDim>               resTotal(cellShape);
      std::vector<par::Array<FluxRes,nDim>>  resStage=par::vec_of_Arrays<FluxRes,nDim>(rk.nstages,cellShape);

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

   // integrate forward in time
      for( size_t tstep=0; tstep<nt; tstep++ )
     {
         Real lmax{};
         for( unsigned int stg=0; stg<rk.nstages; stg++ )
        {
         // calculate differences
            dq=gradientCalcP( mesh.cells, q1 );

         // accumulate flux residual
            resStage[stg] = residualCalc( mesh, hoflux, q1, dq );

         // calculate maximum stable timestep for this timestep
            if( stg==0 ){ lmax = spectralRadius( mesh.cells, resStage[stg] ); }

         // accumulate stage residual
            resTotal = rungeKuttaAccumulation( rk, stg, resStage );

         // integrate cell residuals forward by dt and average over cell volume
            q1 = eulerForwardUpdate( mesh.cells, species, rk.beta[stg]*cfl, lmax, resTotal, q );
        }
         par::copy( q, q1 );
     }

   // write solution to file
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
               writeState( solutionFile, mesh.cells[{i,j}], species, set2State( species, q[{i,j}] ) );
           }
            solutionFile << std::endl;
        }

/*
         par::for_each( // write state to file
                        [&]( const SolVarSet& q0 ) -> void
                           { writeState( solutionFile, species, set2State( species, q0 ) ); },
                        // solution array
                        q );
*/
     }
      else
     {
         std::cout << "cannot open \"data/twoD/periodic/scalarAdvection/result.dat\" for writing\n" << std::endl;
         return 1;
     }
      solutionFile.close();

      return 0;
  }

