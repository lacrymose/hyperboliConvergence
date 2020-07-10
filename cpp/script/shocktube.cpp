
# include <ode.h>

# include <spatial/gradientCalc.h>
# include <spatial/residualCalc.h>
# include <spatial/spectralRadius.h>
# include <spatial/eulerForwardUpdate.h>

# include <limiters/limiter.h>

# include <conservationLaws/euler/euler.h>
# include <conservationLaws/euler/boundaryConditions.h>

# include <geometry/geometry.h>

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

constexpr ShockTube1D problem = ShockTube1D::LowMach;

using Real = double;

constexpr int  nx  = 128;
constexpr int  nt  = 192;
constexpr Real cfl = 0.80;

constexpr Real dt = 0.5e-3;

constexpr BasisT SolutionBasis = BasisT::Primitive;

//using Flux = RusanovFlux<Law>;

//using Flux = AusmPlusUP<LowMachScaling::Acoustic,
//                        LowMachScaling::Acoustic>;

using Flux = Slau<LowMachScaling::Acoustic,
                  LowMachScaling::Acoustic>;

using Limiter = Limiters::VanAlbada2;


// ------ typedefs -------------------

using SolVarSet     = VariableSet<  Law,nDim,SolutionBasis, Real>;
using SolVarDelta   = VariableDelta<Law,nDim,SolutionBasis, Real>;

using FluxRes = FluxResult<Law,nDim,Real>;

using Point   = geom::Point<  nDim,Real>;
using Surface = geom::Surface<nDim,Real>;
using Volume  = geom::Volume< nDim,Real>;


// ------- i/o ------- 

   void writeState( std::ofstream& os, const Species<Law,Real> species, const State<Law,1,Real>& state )
  {
      const Real gamma = species.gamma;
      const Real a = sqrt( state.speedOfSound2() );
      const Real r = state.density();
      const Real u = state.velocity(0);
      const Real p = state.pressure();
      const Real s = a*a / ( gamma*pow(r,gamma-1) );

      os << r << " "
         << u/a << " "
         << p << " "
         << s << std::endl;
  }


// ------- run script ------- 

   int main()
  {
      const par::Shape<nDim> cellShape{nx};

      const ODE::Explicit::RungeKutta<Real> rk = ODE::Explicit::ssp33<Real>();

      std::cout << std::scientific;
      std::cout.precision(8);

   // setup
      const Species<Law,Real> species = get_air_species<Real>();
      const Flux    flux{};
      const Limiter limiter{};

   // initialise mesh
      const geom::Mesh<nDim,Real> mesh = geom::make_linspace_mesh<Real>( cellShape, 0,nx );

   // initialise solution to left/right states
      par::Array<SolVarSet,nDim> q  = shocktube_initial_solution<SolVarSet>( problem, species, mesh.cells );
      par::Array<SolVarSet,nDim> q1 = par::copy( q );
      par::Array<SolVarSet,nDim> q2 = par::copy( q );

   // gradient array
      par::Array<SolVarDelta,nDim>   dq(cellShape);

   // residual evaluated and accumulated at each runge-kutta stage
      // MDArray has no copy constructor so must construct in-place
      par::Array<FluxRes,nDim>               resTotal(cellShape);
      std::vector<par::Array<FluxRes,nDim>>  resStage;
      for( int i=0; i<rk.nstages; i++ ){ resStage.emplace_back(cellShape); }

   // boundary values
      const std::array<SolVarSet,2> qb0 = shocktube_initial_leftright<SolVarSet>( problem, species );

      std::cout << "left  state: " << qb0[0] << std::endl;
      std::cout << "right state: " << qb0[1] << std::endl;

   // high order reconstruction and flux functions
      const auto hoflux = [&species, &limiter, &flux]
                          ( const Surface&      face,
                            const SolVarDelta& gradl,
                            const SolVarDelta& gradr,
                            const SolVarSet&      ql,
                            const SolVarSet&      qr ) -> FluxRes
     {
      // central and biased differences
         const SolVarDelta dqc = qr - ql;
         const SolVarDelta dql = gradl - dqc;
         const SolVarDelta dqr = gradr - dqc;

      // limited differences
         const SolVarDelta Slopel = limiter( dqc, dql );
         const SolVarDelta Sloper = limiter( dqc, dqr );

      // inviscid flux
         return flux( species, face, ql+0.5*Slopel,
                                     qr-0.5*Sloper );
     };

   // integrate forward in time
      for( int tstep=0; tstep<nt; tstep++ )
     {
         Real lmax{};
         for( int j=0; j<rk.nstages; j++ )
        {
         // update boundary conditions using riemann invariants
            const std::array<SolVarSet,2> qbc = updateBCs( species, mesh.nodes, q1, qb0 );

         // calculate differences
            gradientCalc( mesh.cells, q1, dq );
   
         // accumulate flux residual
            residualCalc( mesh, hoflux, qbc, q1, dq, resStage[j] );
   
         // calculate maximum stable timestep for this timestep
            if( j==0 )
           {
               lmax = spectralRadius( mesh.cells, resStage[j] );
               lmax = cfl/dt;
//             std::cout << cfl/lmax << std::endl;
           }

         // accumulate stage residual
            par::fill( resTotal, FluxRes{} );

            for( size_t i=0; i<nx; i++ )
           {
               for( int k=0; k<=j; k++ )
              {
                  resTotal[{i}].flux+=rk.alpha[j][k]*resStage[k][{i}].flux;
              }
           }
   
         // integrate cell residuals forward by dt and average over cell volume
            eulerForwardUpdate( mesh.cells, species, rk.beta[j]*cfl, lmax, resTotal, q, q2 );
            std::swap( q1,q2 );
        }
         par::copy( q, q1 );
     }

   // write solution to file
      std::ofstream solutionFile = []()
     {
         if( problem == ShockTube1D::Sods )
        {
            return std::ofstream( "data/shocktube/sods/result.dat" );
        }
         else if( problem == ShockTube1D::LowMach )
        {
            return std::ofstream( "data/shocktube/lowmach/result.dat" );
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
                        q );
     }
      else
     {
         std::cout << "cannot open \"data/riemannProblem/result.dat\" for writing\n" << std::endl;
         return 1;
     }
      solutionFile.close();

      return 0;
  }

