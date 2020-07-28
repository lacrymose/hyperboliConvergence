
# include <ode.h>

# include <spatial/gradientCalc.h>
# include <spatial/residualCalc.h>
# include <spatial/spectralRadius.h>
# include <spatial/eulerForwardUpdate.h>

# include <limiters/limiter.h>
# include <conservationLaws/euler/euler.h>
# include <conservationLaws/euler/boundaryConditions.h>

# include <mesh/mesh.h>
# include <mesh/generate/oneD.h>

# include <parallalg/algorithm.h>
# include <parallalg/array.h>

# include <vector>
# include <iostream>
# include <fstream>

# include <cases/euler/oneD/periodic.h>

/*
 * One dimensional periodic test cases for the Euler equations
 */

constexpr LawType Law = LawType::Euler;
constexpr int nDim = 1;
using BasisT = BasisType<Law>;
using Real = long double;

// ------- User Inputs -------

constexpr Periodic1D problem = Periodic1D::Soundwave;
//constexpr Periodic1D problem = Periodic1D::AcousticEntropy;

constexpr Real mach = 1.e-8;

// isolated soundwave
constexpr Real amplitude = 0.1;
constexpr unsigned int nwavelengths=1;

// acoustic - entropy interaction
constexpr Real centre_a = 0.2;
constexpr Real centre_s = 0.6;

constexpr Real width_a = 0.02;
constexpr Real width_s = 0.04;

constexpr Real amplitude_a = 0.05;
constexpr Real amplitude_s = 0.00*mach;

// discretisation
constexpr int  nx  = 32;
constexpr int  nt  = 128*128+40;
constexpr Real cfl = 0.5;

constexpr BasisT SolutionBasis = BasisT::Primitive;
constexpr BasisT   ExtrapBasis = BasisT::Primitive;

//using Flux = CentralFlux<Law>;
//using Flux = RoeFlux<Law>;

//using Flux = AusmPlusUP<LowMachScaling::Convective,
//                        LowMachScaling::Acoustic>;

using Flux = Slau<LowMachScaling::Convective,
                  LowMachScaling::Acoustic>;

using Limiter = Limiters::MonotonizedCentral2;
//using Limiter = Limiters::NoLimit3;


// ------ typedefs -------------------

using SolVarSet     = VariableSet<  Law,nDim,SolutionBasis, Real>;
using SolVarDelta   = VariableDelta<Law,nDim,SolutionBasis, Real>;

using ExtrapDelta   = VariableDelta<Law,nDim,ExtrapBasis, Real>;

using FluxRes = FluxResult<Law,nDim,Real>;

using MeshT = Mesh<nDim,Real>;
using Face  = MeshT::Face;


// ------- i/o -------

   void writeState( std::ofstream& os,
                    const Species<Law,  Real>& species,
                    const State<  Law,1,Real>&   state )
  {
      const Real gam = species.gamma;

      const Real r0 = 1.;
      const Real u0 = 1.;
      const Real p0 = r0/(mach*mach*gam);
      const Real a0 = u0/mach;
      const Real s0 = a0*a0 / ( gam*pow(r0,gam-1.) );

      const Real a = sqrt( state.speedOfSound2() );
      const Real r = state.density();
      const Real u = state.velocity(0);
      const Real p = state.pressure();
      const Real s = a*a / ( gam*pow(r,gam-1.) );

      os << (r-r0)/r0 << " "
         << (u-u0)/u0 << " "
         << (p-p0)/p0 << " "
         << (s-s0)/s0 << std::endl;
  }


// ------- run script -------

   int main()
  {
      const par::Shape<nDim> cellShape{nx};

      const ODE::Explicit::RungeKutta<Real> rk = ODE::Explicit::ssp33<Real>();

      std::cout << std::scientific;
      std::cout.precision(8);

   // setup
      const Species<Law,Real> species = []()
     {
         Species<Law,Real> gas = get_air_species<Real>();
         gas.minf=mach;
         return gas;
     }();
      const Flux    flux{};
      const Limiter limiter{};

   // initialise mesh
      const MeshT mesh = make_linspace_mesh<Real>( cellShape, 0,nx );

   // initialise solution
      par::Array<SolVarSet,nDim> q = [&]() -> par::Array<SolVarSet,nDim>
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

      par::Array<SolVarSet,nDim> q1 = par::copy( q );
      par::Array<SolVarSet,nDim> q2 = par::copy( q );

   // gradient array
      par::Array<SolVarDelta,nDim>   dq(cellShape);

   // residual evaluated and accumulated at each runge-kutta stage
      // MDArray has no copy constructor so must construct in-place
      par::Array<FluxRes,nDim>               resTotal(cellShape);
      std::vector<par::Array<FluxRes,nDim>>  resStage;
      for( int i=0; i<rk.nstages; i++ ){ resStage.emplace_back(cellShape); }

   // high order reconstruction and flux functions
      const auto hoflux = [&species, &limiter, &flux]
                          ( const Face&          face,
                            const SolVarDelta&  gradl,
                            const SolVarDelta&  gradr,
                            const SolVarSet&       ql,
                            const SolVarSet&       qr ) -> FluxRes
     {
      // average state
         const State<Law,nDim,Real> avg = set2State( species, ql+0.5*(qr-ql) );

      // central and l/r biased differences in solution basis
         const SolVarDelta dqc = qr - ql;
         const SolVarDelta dql = gradl - dqc;
         const SolVarDelta dqr = gradr - dqc;

      // central and l/r biased differences in extrapolation basis aligned with face
         const ExtrapDelta dvc = delta2Delta<ExtrapDelta>( species, avg, rotateToMetric( face.metric, dqc ) );
         const ExtrapDelta dvl = delta2Delta<ExtrapDelta>( species, avg, rotateToMetric( face.metric, dql ) );
         const ExtrapDelta dvr = delta2Delta<ExtrapDelta>( species, avg, rotateToMetric( face.metric, dqr ) );

      // limited differences in solution basis aligned with background coordinate system
         const SolVarDelta slopel = rotateFromMetric( face.metric,
                                                      delta2Delta<SolVarDelta>( species,
                                                                                avg,
                                                                                limiter( dvc, dvl ) ) );
         const SolVarDelta sloper = rotateFromMetric( face.metric,
                                                      delta2Delta<SolVarDelta>( species,
                                                                                avg,
                                                                                limiter( dvc, dvr ) ) );
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

