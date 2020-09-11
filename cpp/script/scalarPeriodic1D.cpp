
# include <ode.h>

# include <spatial/gradientCalc.h>
# include <spatial/residualCalc.h>
# include <spatial/spectralRadius.h>
# include <spatial/eulerForwardUpdate.h>

# include <limiters/limiter.h>
# include <conservationLaws/scalarAdvection/scalarAdvection.h>
# include <conservationLaws/scalarAdvection/boundaryConditions.h>

# include <solutionField/solutionField.h>

# include <mesh/mesh.h>
# include <mesh/generate/oneD.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <vector>
# include <iostream>
# include <fstream>

# include <cmath>

/*
 * One dimensional periodic test cases for the Euler equations
 */

constexpr LawType Law = LawType::ScalarAdvection;
constexpr int nDim = 1;
using BasisT = BasisType<Law>;
using Real = long double;

// ------- User Inputs -------

constexpr Periodic1D problem = Periodic1D::Soundwave;

// monochromatic soundwave
constexpr Real amplitude = 0.1;
constexpr unsigned int nwavelengths=1;
constexpr Real velocity = 1.;

// discretisation
constexpr int  nx  = 32;
constexpr int  nt  = 128;
constexpr Real cfl = 0.5;

constexpr BasisT SolutionBasis = BasisT::Conserved;

using Flux = RusanovFlux<Law>;

//using Limiter = Limiters::MonotonizedCentral2;
using Limiter = Limiters::NoLimit3;


// ------ typedefs -------------------

using SolVarSet   = VariableSet<  Law,nDim,SolutionBasis,Real>;
using SolVarDelta = VariableDelta<Law,nDim,SolutionBasis,Real>;
using SolField    = SolutionField<Law,nDim,SolutionBasis,Real>;

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

      std::cout << std::scientific;
      std::cout.precision(8);

   // setup
      const Species<Law,Real> species{};
      const Flux                 flux{};
      const Limiter           limiter{};

   // initialise mesh
      const MeshT mesh = make_linspace_mesh<Real>( cellShape, 0,nx );

   // create solution field
      SolField q( cellShape);

   // initialise interior to sinusoid
      const auto initial_conditions = [&]( const Mesh::Cell& c ) -> SolVarSet
     {
         const Real x = c.centre/nx;
         const Real l = nwavelengths*x;
         return { velocity, amplitude*sin(2.*M_PI*l) };
     };
      par::transform( initial_conditions, q.interior, mesh.cells );

   // initialise periodic boundaries
      for( SolField::BCType&  bc : q.bcTypes ){ bc = SolField::BCType::Periodic; }
      for( SolField::VarField& v : q.boundaries ) { v[{0}]=0.; v[{1}]=0.; }

      SolField q1 = copy(q);
      SolField q2 = copy(q);

   // gradient array
      par::Array<SolVarDelta,nDim> dq(cellShape);

   // residual evaluated and accumulated at each runge-kutta stage
      par::Array<FluxRes,nDim> resTotal(cellShape);
      auto resStage = par::vec_of_Arrays<FluxRes,nDim>( rk.nstages, cellShape );

   // high order reconstruction and flux functions
      const auto hoflux = [&species, &limiter, &flux]
                          ( const Face&          face,
                            const SolVarDelta&  gradl,
                            const SolVarDelta&  gradr,
                            const SolVarSet&       ql,
                            const SolVarSet&       qr ) -> FluxRes
     {
      // central and l/r biased differences in solution basis with background coordinates
         const SolVarDelta dqc = qr - ql;
         const SolVarDelta dql = gradl - dqc;
         const SolVarDelta dqr = gradr - dqc;

      // limited differences
         const SolVarDelta slopel = limiter( dqc, dql );
         const SolVarDelta sloper = limiter( dqc, dqr );

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
            gradientCalc( mesh.cells, q1, dq );

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
         copy( q, q1 );
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

