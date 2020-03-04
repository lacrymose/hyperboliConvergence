
# include <ode.h>
# include <controls.h>
# include <array1D/array1D.h>
# include <idealGas2D/idealGas2D.h>
# include <fieldOperations/fieldOperations.h>

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
      Array::Array1D<  SolutionVariables>  qm;
      Array::Array1D<  SolutionVariables>  q[3];
      Array::Array1D<Gas::ConservedDelta>  r[2];
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

      bool contact=false;

      int      i,j,k;

      float  lmax,dt,dtau;

      std::cout << std::scientific;
      std::cout.precision(2);

      gas.air();
//    outerTimeDerivative.eulerBackward1();
//    outerTimeDerivative.backwardDifference2();
      outerTimeDerivative.trapeziumRule2();
      explicitIntegrator.ssp34();

   // read parameters
      if( initialStatesFile.is_open() )
     {
         for( i=0; i<4; i++ ){ initialStatesFile >> qvl[i]; }
         for( i=0; i<4; i++ ){ initialStatesFile >> qvr[i]; }
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

      q[0].resize(grid.n);
      q[1].resize(grid.n);
      q[2].resize(grid.n);
      qm.resize(grid.n);

      r[0].resize(grid.n);
      r[1].resize(grid.n);
      s.resize( grid.n);

      contact=true;
      if( contact )
     {
         float mach=0.2;
         float perturbation=0.05;

         float temp,press;

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

      for( i=0;        i<grid.n/4; i++ ){ q[0][i] = ql; }
      for( i=grid.n/4; i<grid.n;   i++ ){ q[0][i] = qr; }

   // timesteps
      q[1]=q[0];
      q[2]=q[0];
      qm  =q[0];
      r[0]=0.;
      r[1]=0.;
      s   =0.;

      for( i=0; i<outerTime.nt; i++ )
     {
      // outer timestep calculation
         lmax=-1;
         fluxResidual( gas, grid.boundaryCondition, flux, q[1],r[1], lmax );
         if( i==0 ){ dt = outerTime.cfl/lmax; }

      // dualtime iterations
         for( j=0; j<innerTime.nt; j++ )
        {
            for( k=0; k<grid.n; k++ )
           {
               s[k] =  outerTimeDerivative.beta[0]*Gas::ConservedDelta( Gas::ConservedVariables( gas, q[0][k] ) )
                     + outerTimeDerivative.beta[1]*Gas::ConservedDelta( Gas::ConservedVariables( gas, q[1][k] ) );
                     + outerTimeDerivative.beta[2]*Gas::ConservedDelta( Gas::ConservedVariables( gas, q[2][k] ) );
               s[k]/= dt;
           }

            lmax=-1;
            fluxResidual( gas, grid.boundaryCondition, flux, q[0],r[0], lmax );
            dtau = innerTime.cfl/lmax;

            for( k=0; k<grid.n; k++ )
           {
               r[0][k] =  outerTimeDerivative.gamma[0]*r[0][k]
                        + outerTimeDerivative.gamma[1]*r[1][k]
                        - s[k];
           }

            dtau/= ( 1. + outerTimeDerivative.beta[0]*dtau/dt );

            eulerForwardLinearUpdate( gas, dtau, q[0],qm, r[0], res );
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

         res=0;
         for( k=0; k<grid.n; k++ )
        {
            for( int m=0; m<4; m++ )
           {
               res[m]+= fabs( Gas::ConservedVariables( gas, q[0][k] )[m]
                             -Gas::ConservedVariables( gas, q[1][k] )[m] );
           }
        }

         std::cout << std::endl;
         std::cout << std::left << std::setw(8) << i;
         std::cout << res << std::endl;
         std::cout << std::endl;

         qm=q[0];
         float relax=0.0;
         for( k=0; k<grid.n; k++ )
        {
            qm[k]= q[0][k] + relax*outerTimeDerivative.beta[0]*SolutionDelta( q[0][k] )
                           + relax*outerTimeDerivative.beta[1]*SolutionDelta( q[1][k] )
                           + relax*outerTimeDerivative.beta[2]*SolutionDelta( q[2][k] );
        }

         q[2]=q[1];
         q[1]=q[0];
         q[0]=qm;
     }

   // save solution
      std::ofstream solutionFile( "data/sods/sods.dat" );
      Gas::State state;
      if( solutionFile.is_open() )
     {
         for( i=0; i<grid.n; i++ )
        {
            state = Gas::State( gas, q[0][i] );
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
