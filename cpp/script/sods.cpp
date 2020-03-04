
# include <ode.h>
# include <controls.h>
# include <limiter.h>
# include <array1D/array1D.h>
# include <idealGas2D/idealGas2D.h>
# include <timestepping/timestepping.h>

# include <iostream>
# include <fstream>
# include <cassert>

namespace Gas = IdealGas2D;

typedef Gas::Conserved        SolutionType;
//typedef Gas::Viscous          SolutionType;

typedef Gas::Rusanov        Flux;
//typedef Gas::Ausm           Flux;
//typedef Gas::Slau           Flux;

typedef Limiters::NoLimit1 Limiter;
//typedef Limiters::MinMod2 Limiter;
//typedef Limiters::Cada3 Limiter;

typedef Gas::VariableSet<  SolutionType>         SolutionVariables;
typedef Gas::VariableDelta<SolutionType>             SolutionDelta;

   int main()
  {
   // solution vector
      Array::Array1D<  SolutionVariables> q;

   // initial conditions
      Gas::ViscousVariables   qvl,qvr;
      SolutionVariables        ql, qr;

      Gas::Species   gas;
      Flux          flux;

      Limiters::VectorFaceLimiter<Limiter>      limiter;

      Controls::GridControls1D                     grid;
      Controls::TimeSteppingControls  outerTimeControls;
      Controls::TimeSteppingControls  innerTimeControls;

      ODE::Implicit::MultiStep       implicitIntegrator;
      ODE::Explicit::RungeKutta      explicitIntegrator;

      bool print=true;
      int  i;

      bool shocktube=true;
      bool contact  =false;

      std::cout << std::scientific;
      std::cout.precision(2);

      gas.air();
//    implicitIntegrator.eulerBackward1();
//    implicitIntegrator.backwardDifference2();
      implicitIntegrator.trapeziumRule2();
      explicitIntegrator.ssp34();

   // read problem parameters
      std::ifstream parametersFile( "data/sods/parameters.dat" );
      if( parametersFile.is_open() )
     {
         parametersFile >>              grid;
         parametersFile >> outerTimeControls;
         parametersFile >> innerTimeControls;
     }
      else
     {
         std::cout << "cannot open \"parameters.dat\" for reading\n" << std::endl;
         return 1;
     }
      parametersFile.close();

   // initialise solution
      q.resize( grid.n );

      if( shocktube )
     {
         grid.boundaryCondition='a';
         gas.minf=1.;
         std::ifstream initialStatesFile( "data/sods/initial.dat" );
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
     }

      if( contact )
     {
         grid.boundaryCondition='p';
         float mach=0.1;
         float perturbation=0.05;

         float temp,press;

         gas.minf = mach;

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

      for( i=0;        i<(int)grid.n/4; i++ ){ q[i] = ql; }
      for( i=grid.n/4; i<(int)grid.n;   i++ ){ q[i] = qr; }


//    TimeStepping::explicitEuler( outerTimeControls,                    grid, gas,flux,limiter, q, print );
//    TimeStepping::rungeKutta(    outerTimeControls,explicitIntegrator, grid, gas,flux,limiter, q, print );
      TimeStepping::dualTime(      outerTimeControls,implicitIntegrator,
                                   innerTimeControls,explicitIntegrator, grid, gas,flux,limiter, q, print );

   // save solution
      std::ofstream solutionFile( "data/sods/sods.dat" );
      Gas::State state;
      if( solutionFile.is_open() )
     {
         for( i=0; i<(int)grid.n; i++ )
        {
            state = Gas::State( gas, q[i] );
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
