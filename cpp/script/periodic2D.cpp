
# include <ode.h>

# include <spatial/gradientCalc.h>
# include <spatial/residualCalc.h>
# include <spatial/spectralRadius.h>
# include <spatial/eulerForwardUpdate.h>

# include <limiters/limiter.h>
# include <conservationLaws/scalarAdvection/scalarAdvection.h>

# include <geometry/geometry.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <vector>
# include <iostream>
# include <fstream>

# include <cases/scalarAdvection/twoD/periodic.h>

/*
 * One dimensional periodic test cases for the Euler equations
 */

constexpr LawType Law = LawType::ScalarAdvection;
constexpr int nDim = 2;
using BasisT = BasisType<Law>;
using Real = float;

// ------- User Inputs -------

// Convecting velocity - angle and magnitude
constexpr Real theta = 45.;
constexpr Real speed = 1.0;

// width of scalar top-hat (cells)
constexpr int nh = 8;

// discretisation
constexpr int  nx  = 24;
constexpr int  nt  = 0;
constexpr Real cfl = 0.5;

constexpr BasisT SolutionBasis = BasisT::Conserved;

using Flux = RusanovFlux<Law>;

using Limiter = Limiters::NoLimit1;


// ------ typedefs -------------------

using SolVarSet  = VariableSet<Law,nDim,SolutionBasis, Real>;
using SolVarDel  = VariableDelta<Law,nDim,SolutionBasis, Real>;
using SolVarGrad = std::array<SolVarDel,nDim>;

using FluxRes = FluxResult<Law,nDim,Real>;

using Point   = geom::Point<  nDim,Real>;
using Surface = geom::Surface<nDim,Real>;
using Volume  = geom::Volume< nDim,Real>;


// ------- i/o -------

   void writeState( std::ofstream& os, const Species<Law,Real> species, const State<Law,nDim,Real>& state )
  {
      const Real u = state.velocity(0);
      const Real v = state.velocity(1);
      const Real q = state.scalar();

      os << u << " "
         << v << " "
         << q << std::endl;
  }


// ------- run script -------

   int main()
  {
      const par::Shape<nDim> cellShape{nx,nx};

      const ODE::Explicit::RungeKutta<Real> rk = ODE::Explicit::ssp33<Real>();

      std::cout << std::scientific;
      std::cout.precision(8);

   // setup
      const Species<Law,Real> species{};
      const Flux                 flux{};
      const Limiter           limiter{};

   // initialise mesh
      const geom::Mesh<nDim,Real> mesh = geom::make_linspace_mesh<Real>( cellShape, 0,nx, 0,nx ); // ???

   // initialise solution
      par::Array<SolVarSet,nDim> q = // ???

      par::Array<SolVarSet,nDim> q1 = par::copy( q );
      par::Array<SolVarSet,nDim> q2 = par::copy( q );

   // gradient array
      par::Array<SolVarGrad,nDim>   dq(cellShape);

   // residual evaluated and accumulated at each runge-kutta stage
      // MDArray has no copy constructor so must construct in-place
      par::Array<FluxRes,nDim>               resTotal(cellShape);
      std::vector<par::Array<FluxRes,nDim>>  resStage;
      for( int i=0; i<rk.nstages; i++ ){ resStage.emplace_back(cellShape); }

   // high order reconstruction and flux functions
      const auto hoflux = [&species, &limiter, &flux]
                          ( const Surface&    face,
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
      for( int tstep=0; tstep<nt; tstep++ )
     {
         Real lmax{};
         for( int j=0; j<rk.nstages; j++ )
        {
         // calculate differences
            gradientCalcP( mesh.cells, q1, dq );

         // accumulate flux residual
            residualCalc( mesh, hoflux, q1, dq, resStage[j] );

         // calculate maximum stable timestep for this timestep
            if( j==0 ){ lmax = spectralRadius( mesh.cells, resStage[j] ); }

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
      std::ofstream solutionFile = []() -> std::ofstream
     {
         if( problem == Periodic1D::Soundwave )
        {
            return std::ofstream( "data/periodic/soundwave/result.dat" );
        }
         else if( problem == Periodic1D::AcousticEntropy )
        {
            return std::ofstream( "data/periodic/acoustic_entropy/result.dat" );
        }
     }();

      if( solutionFile.is_open() )
     {
         solutionFile << std::scientific;
         solutionFile.precision(12);

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

