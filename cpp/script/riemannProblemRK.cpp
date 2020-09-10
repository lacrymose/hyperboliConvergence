
# include <ode.h>

# include <spatial/gradientCalc.h>
# include <spatial/residualCalc.h>
# include <spatial/spectralRadius.h>
# include <spatial/eulerForwardUpdate.h>

# include <limiters/limiter.h>

# include <conservationLaws/scalarAdvection/scalarAdvection.h>
# include <conservationLaws/euler/euler.h>
# include <conservationLaws/euler/boundaryConditions.h>

# include <geometry/geometry.h>

# include <mdarray/mdarray.h>

# include <vector>
# include <iostream>
# include <fstream>

# include <cassert>


// ------- Inputs ------- 
using Real = float;

constexpr LawType Law = LawType::Euler;
//constexpr LawType Law = LawType::ScalarAdvection;
constexpr int nDim = 1;

constexpr int  nx  = 128;
constexpr int  nt  = 192;
constexpr Real cfl = 1.60;

using BasisT = BasisType<Law>;
//constexpr BasisT SolutionBasis = BasisT::Conserved;
constexpr BasisT SolutionBasis = BasisT::Viscous;

using ConvectiveFlux  = CentralFlux<Law>;
using DissipativeFlux = RusanovDissipation<Law>;

using ConvectiveLimiter  = Limiters::NoLimit3;
using DissipativeLimiter = Limiters::Cada3;

constexpr bool print=false;

// -----------------------------------

constexpr BasisT ConservedBasis = BasisT::Conserved;

using SolutionVarSet    = VariableSet<  Law,nDim,SolutionBasis, Real>;
using SolutionVarDelta  = VariableDelta<Law,nDim,SolutionBasis, Real>;
using ConservedVarDelta = VariableDelta<Law,nDim,ConservedBasis,Real>;

using FluxRes = FluxResult<Law,nDim,Real>;

using Point   = geom::Point<  nDim,Real>;
using Surface = geom::Surface<nDim,Real>;
using Volume  = geom::Volume< nDim,Real>;


# include "../script/riemannProblem.h"

   int main()
  {
      const Dims<nDim> cellDims{nx};
      const Dims<nDim> nodeDims{nx+1};

      const ODE::Explicit::RungeKutta<Real> rk = ODE::Explicit::ssp34<Real>();

      std::cout << std::scientific;
      std::cout.precision(3);

   // setup
      const Species<Law,Real> species = get_species( ConservedBasis );
      const ConvectiveFlux    cFlux{};
      const ConvectiveLimiter cLimiter{};

      const DissipativeFlux    dFlux{};
      const DissipativeLimiter dLimiter{};

   // solutions arrays
      MDArray<SolutionVarSet,  nDim>    q(cellDims);
      MDArray<SolutionVarSet,  nDim>   q1(cellDims);
      MDArray<SolutionVarSet,  nDim>   q2(cellDims);
      MDArray<SolutionVarDelta,nDim>   dq(cellDims);

   // residual evaluated and accumulated at each runge-kutta stage
      // MDArray has no copy constructor so must construct in-place
      MDArray<FluxRes,nDim>               resTotal(cellDims);
      std::vector<MDArray<FluxRes,nDim>>  resStage;
      for( int i=0; i<rk.nstages; i++ ){ resStage.emplace_back(cellDims); }

   // initialise mesh
      const MDArray<Point,nDim> nodes = [nodeDims]()
     {
         MDArray<Point,nDim> nds(nodeDims);
         int i=0;
         for( Point& p : nds.elems ){ p = Point{ Real(i++) }; };
         return nds;
     }();
      const MDArray<Volume,1> cells = geom::dual( nodes );

   // initialise solution to left/right states
      const SolutionVarSet ql0 = set2Set<SolutionVarSet>( species, initialLeft(  SolutionBasis ) );
      const SolutionVarSet qr0 = set2Set<SolutionVarSet>( species, initialRight( SolutionBasis ) );

      std::cout << "left  state: " << ql0 << std::endl;
      std::cout << "right state: " << qr0 << std::endl;

      for( size_t i=0;    i<nx/2; i++ ){ q[{i}]=ql0; q1[{i}]=ql0; q2[{i}]=ql0; }
      for( size_t i=nx/2; i<nx;   i++ ){ q[{i}]=qr0; q1[{i}]=qr0; q2[{i}]=qr0; }

   // boundary values
      std::array<SolutionVarSet,2> qb0{ql0,qr0};

   // high order reconstruction and flux functions
      const auto hoflux = [&species, &cLimiter, &cFlux, &dLimiter, &dFlux]
                          ( const Surface&          face,
                            const SolutionVarDelta&  dql,
                            const SolutionVarDelta&  dqc,
                            const SolutionVarDelta&  dqr,
                            const SolutionVarSet&     ql,
                            const SolutionVarSet&     qr ) -> FluxRes
     {
      // convective flux
         const SolutionVarDelta cSlopel = cLimiter( dqc, dql );
         const SolutionVarDelta cSloper = cLimiter( dqc, dqr );

         const FluxRes cfr = cFlux( species, face, ql+0.5*cSlopel,
                                                   qr-0.5*cSloper );

      // dissipative flux
         const SolutionVarDelta dSlopel = dLimiter( dqc, dql );
         const SolutionVarDelta dSloper = dLimiter( dqc, dqr );

         const FluxRes dfr = dFlux( species, face, ql+0.5*dSlopel,
                                                   qr-0.5*dSloper );

      // accumulated flux
         return FluxRes( cfr.flux + dfr.flux,
                         std::max( cfr.lambda, dfr.lambda ) );
     };

   // integrate forward in time
      for( int t=0; t<nt; t++ )
     {
         Real lmax{};
         for( int j=0; j<rk.nstages; j++ )
        {
         // update boundary conditions using riemann invariants
            const std::array<SolutionVarSet,2> qbc = updateBCs( species, nodes, q1, qb0 );

         // calculate differences
            gradientCalc( cells, q1, dq );
   
         // accumulate flux residual
            residualCalc( cells, nodes, hoflux, qbc, q1, dq, resStage[j] );
   
         // calculate maximum stable timestep for this timestep
            if( j==0 ){ lmax = spectralRadius( cells, resStage[j] ); }

         // accumulate stage residual
            for( FluxRes& r : resTotal.elems ){ r = FluxRes{}; }
            for( size_t i=0; i<nx; i++ )
           {
               for( int k=0; k<=j; k++ )
              {
                  resTotal[{i}].flux+=rk.alpha[j][k]*resStage[k][{i}].flux;
              }
           }
   
         // integrate cell residuals forward by dt and average over cell volume
            eulerForwardUpdate( cells, species, rk.beta[j]*cfl, lmax, resTotal, q, q2 );
            std::swap( q1,q2 );
        }
         for( size_t i=0; i<nx; i++ ){ q[{i}]=q1[{i}]; }
     }

      std::ofstream solutionFile( "data/riemannProblem/result.dat" );
      if( solutionFile.is_open() )
     {
         for( const SolutionVarSet& q0 : q.elems )
        {
            writeState( solutionFile, set2State( species, q0 ) );
        }
     }
      else
     {
         std::cout << "cannot open \"data/riemannProblem/result.dat\" for writing\n" << std::endl;
         return 1;
     }
      solutionFile.close();

      return 0;
  }

