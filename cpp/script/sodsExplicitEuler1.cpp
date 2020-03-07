
# include <array1D/array1D.h>
# include <idealGas2D/idealGas2D.h>
# include <fieldOperations/fieldOperations.h>


# include <iostream>
# include <fstream>
# include <cassert>
# include <cmath>


namespace Gas = IdealGas2D;

//typedef Gas::Conserved    SolutionType;
typedef Gas::Viscous      SolutionType;

//typedef Gas::LaxFriedrichs   Flux;
typedef Gas::Ausm            Flux;

typedef Gas::VariableSet<  SolutionType>   SolutionVariables;
typedef Gas::VariableDelta<SolutionType>   SolutionDelta;

   int main()
  {
      std::ifstream initialStatesFile( "data/sods/initial.dat" );
      std::ifstream parametersFile( "data/sods/parameters.dat" );

   // initial conditions
      IdealGas2D::ViscousVariables    qvl,qvr;
      SolutionVariables               ql,qr;

   // flux / conserved state;
      IdealGas2D::ConservedDelta    f,res;

   // initialise arrays
      Array::Array1D<     SolutionVariables> q0;
      Array::Array1D<     SolutionVariables> q1;
      Array::Array1D<Gas::ConservedDelta> r;


      Gas::Species  gas;
      Flux         flux;

      bool contact=false;

      int      nc,nt;
      int          i;

      float      cfl;
      float  lmax,dt;

      gas.air();

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
         parametersFile >>  nc;
         parametersFile >>  nt;
         parametersFile >> cfl;
     }
      else
     {
         std::cout << "cannot open \"parameters.dat\" for reading\n" << std::endl;
         return 1;
     }
      parametersFile.close();

   // ensure domain is symmetric and waves will not reach boundaries
      assert( nc%2==0 );
//    assert( nt*cfl<0.5*nc );

      q0.resize(nc);
      q1.resize(nc);
      r.resize(nc);

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
      for( i=0;    i<nc/4; i++ ){ q0[i] = ql; }
      for( i=nc/4; i<nc;   i++ ){ q0[i] = qr; }

   // timesteps
      r=0.;
      q1=0.;
      for( i=0; i<nt; i++ )
     {
         lmax=-1;

      // flux residual
         fluxResidual( gas, 'a', flux, q0,r, lmax );

         dt = cfl/lmax;

      // update
         eulerForwardNonlinearUpdate( gas, dt, q0,q1, r, res );
         q0=q1;
         std::cout << res << std::endl;
     }

   // save solution
      std::ofstream solutionFile( "data/sods/sods.dat" );
      IdealGas2D::State s;
      if( solutionFile.is_open() )
     {
         for( i=0; i<nc; i++ )
        {
            s = IdealGas2D::State( gas, q0[i] );
            solutionFile << s.density()     << " "
                         << s.velocityX()   << " "
                         << s.velocityY()   << " "
                         << s.pressure()    << " "
                         << s.temperature() << std::endl;
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
