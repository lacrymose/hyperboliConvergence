
# include <ode.h>
# include <controls.h>
# include <limiter.h>
# include <array1D/array1D.h>
# include <idealGas2D/idealGas2D.h>
# include <fieldOperations/fieldOperations.h>

# include <types.h>

# include <iostream>
# include <fstream>
# include <assert.h>
# include <cmath>


namespace Gas = IdealGas2D;

typedef Gas::Conserved    SolutionType;
//typedef Gas::Viscous      SolutionType;

typedef Gas::Rusanov   Flux;
//typedef Gas::Ausm      Flux;
//typedef Gas::Slau      Flux;

typedef Limiters::NoLimit1 Limiter;
//typedef Limiters::MinMod2 Limiter;
//typedef Limiters::Cada3 Limiter;

typedef Gas::VariableSet<  SolutionType>   SolutionVariables;
typedef Gas::VariableDelta<SolutionType>   SolutionDelta;

   int main()
  {
      std::ifstream initialStatesFile( "data/sods/initial.dat" );
      std::ifstream parametersFile( "data/sods/parameters.dat" );

   // initial conditions
      Gas::ViscousVariables    qvl,qvr;
      SolutionVariables         ql, qr;

   // initialise arrays
      Array::Array1D<  SolutionDelta>      dq;
      Array::Array1D<  SolutionVariables>  qm;
      std::vector<Array::Array1D<  SolutionVariables>>  q;
      std::vector<Array::Array1D<Gas::ConservedDelta>>  r;
      Array::Array1D<Gas::ConservedDelta>  rm;
      Array::Array1D<Gas::ConservedDelta>  s;

      Gas::ConservedDelta    res;
      Gas::ConservedDelta    res0;

      Gas::Species  gas;
      Flux         flux;

      Limiters::VectorFaceLimiter<Limiter>      limiter;

      Controls::GridControls1D                     grid;
      Controls::TimeSteppingControls          outerTime;
      Controls::TimeSteppingControls          innerTime;

      ODE::Implicit::MultiStep       implicitIntegrator;
      ODE::Explicit::RungeKutta      explicitIntegrator;

      bool contact=true;

      size_t    ii,k;
      int        i,j;

      Types::Real  lmax,dt,dtau;

      std::cout << std::scientific;
      std::cout.precision(2);

      gas.air();
//    implicitIntegrator.eulerBackward1();
//    implicitIntegrator.backwardDifference2();
      implicitIntegrator.trapeziumRule2();
      explicitIntegrator.ssp34();

   // read parameters
      if( initialStatesFile.is_open() )
     {
         for( j=0; j<4; j++ ){ initialStatesFile >> qvl[j]; }
         for( j=0; j<4; j++ ){ initialStatesFile >> qvr[j]; }
     }
      else
     {
         std::cout << "cannot open \"initial.dat\" for reading\n" << std::endl;
         return 1;
     }
      initialStatesFile.close();

      if( parametersFile.is_open() )
     {
         parametersFile >>      grid;
         parametersFile >> outerTime;
         parametersFile >> innerTime;
     }
      else
     {
         std::cout << "cannot open \"parameters.dat\" for reading\n" << std::endl;
         return 1;
     }
      parametersFile.close();

   // ensure domain is symmetric and waves will not reach boundaries
      assert( grid.n%2==0 );

      q.resize(implicitIntegrator.nsteps);
      r.resize(implicitIntegrator.nresid);

      for( ii=0; ii<q.size(); ii++ ){ q[ii].resize( grid.n ); }
      for( ii=0; ii<r.size(); ii++ ){ r[ii].resize( grid.n ); }

      rm.resize( grid.n   );
      qm.resize( grid.n   );
      dq.resize( grid.n+2 );

      s.resize( grid.n);

      if( contact )
     {
         float mach=0.2;
         float perturbation=0.05;

         float temp,press;

         grid.boundaryCondition='p';

         press= 1./(gas.gamma*mach*mach);
         temp = press/gas.Rgas;

         qvl[0]=1.;      qvr[0]=1.;
         qvl[1]=0.;      qvr[1]=0.;
         qvl[2]=temp;    qvr[2]=temp;
         qvl[3]=press;   qvr[3]=press;

         qvl[2]*=(1.+perturbation);
     }

      ql = SolutionVariables( gas, qvl );
      qr = SolutionVariables( gas, qvr );

      for( ii=0;        ii<grid.n/4; ii++ ){ q[0][ii] = ql; }
      for( ii=grid.n/4; ii<grid.n;   ii++ ){ q[0][ii] = qr; }

   // timesteps
      for( ii=1; ii<q.size(); ii++ ){ q[ii]=q[0]; }
      for( ii=1; ii<r.size(); ii++ ){ r[ii]=  0.; }
      qm  =q[0];
      s   =0.;

      for( i=0; i<outerTime.nt; i++ )
     {
      // outer timestep calculation
         lmax=-1;
         gradientCalculation( grid, q[1],dq );
         if( r.size()==2 )
        {
            fluxResidual( gas, grid.boundaryCondition, flux,limiter, q[0],dq,r[1], lmax );
        }
         else
        {
            fluxResidual( gas, grid.boundaryCondition, flux,limiter, q[0],dq,r[0], lmax );
        }
         if( i==0 ){ dt = outerTime.cfl/lmax; }

      // dualtime iterations
         for( j=0; j<innerTime.nt; j++ )
        {
            lmax=-1;
            gradientCalculation( grid, q[0],dq );
            fluxResidual( gas, grid.boundaryCondition, flux,limiter, q[0],dq,r[0], lmax );

            dtau = innerTime.cfl/lmax;
            dtau/= ( 1. + implicitIntegrator.beta[0]*dtau/dt );

            dualTimeResidual( implicitIntegrator, gas, grid,
                              dt,q,r, rm );


            eulerForwardLinearUpdate( gas, dtau, q[0],qm, rm, res );
            q[0]=qm;

            if( j==0 )
           {
               res0=res;
               std::cout << std::left << std::setw(4) << i;
               std::cout << std::left << std::setw(4) << j;
               std::cout << res << std::endl;
           }
            else
           {
               for( int m=0; m<4; m++ ){ res[m]/=res0[m]; }
               res[2]=0;
               std::cout << std::left << std::setw(4) << i;
               std::cout << std::left << std::setw(4) << j;
               std::cout << res << std::endl;
           }
        }

         residual( gas, q[0],q[1], res );

         std::cout << std::endl;
         std::cout << std::left << std::setw(8) << i;
         std::cout << res << std::endl;
         std::cout << std::endl;

         for( ii=q.size()-1; ii>0; ii-- ){ q[ii]=q[ii-1]; }
     }

   // save solution
      std::ofstream solutionFile( "data/sods/sods.dat" );
      Gas::State state;
      if( solutionFile.is_open() )
     {
         for( k=0; k<grid.n; k++ )
        {
            state = Gas::State( gas, q[0][k] );
            solutionFile << state.density()     << " "
                         << state.velocityX()   << " "
                         << state.velocityY()   << " "
                         << state.pressure()    << " "
                         << state.temperature() << std::endl;
        }
     }
      else
     {
         std::cout << "cannot open \"sods.dat for writing\"\n" << std::endl;
         return 1;
     }
      solutionFile.close();

      return 0;
  }
