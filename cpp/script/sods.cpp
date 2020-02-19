# include <idealGas2D/idealGas2D.h>

# include <iostream>
# include <fstream>
# include <assert.h>
# include <math.h>

//auto& flux = IdealGas2D::laxFriedrichs;
auto& flux = IdealGas2D::ausm;

typedef IdealGas2D::ConservedVariables SolutionVariables;
//typedef IdealGas2D::ViscousVariables SolutionVariables;

   int main()
  {
      std::ifstream initialStatesFile( "initial.dat" );
      std::ifstream parametersFile( "parameters.dat" );

   // initial conditions
      IdealGas2D::ViscousVariables    qvl,qvr;
      SolutionVariables               ql,qr;

   // flux / conserved state;
      IdealGas2D::ConservedVariables    f,qc;

   // solution and residual vectors
      SolutionVariables                *q;
      IdealGas2D::ConservedVariables   *r;

      IdealGas2D::Species  gas;
      IdealGas2D::State    sl,sr;

      int      nc,nt;
      int        i,j;

      float          cfl;
      float    l,lmax,dt;
      float n[3]={1,0,1};

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
      assert( nt*cfl<0.5*nc );

   // initialise
      q = new SolutionVariables[nc];
      r = new IdealGas2D::ConservedVariables[nc];

      ql = SolutionVariables( gas, qvl );
      qr = SolutionVariables( gas, qvr );
      for( i=0;    i<nc/2; i++ ){ q[i] = ql; }
      for( i=nc/2; i<nc;   i++ ){ q[i] = qr; }

   // timesteps
      for( i=0; i<nt; i++ )
     {
         lmax=-1;

      // flux residual
         for( j=0; j<nc-1; j++ )
        {
            sl = IdealGas2D::State( gas, q[j]   );
            sr = IdealGas2D::State( gas, q[j+1] );

//          flux( gas, n, sl,sr, f,l );
            flux( gas, n,
                  IdealGas2D::State( gas, q[j]   ),
                  IdealGas2D::State( gas, q[j+1] ),
                  f,l );

            lmax = fmax( lmax, l );

            r[j]  -= f;
            r[j+1]+= f;
        }
         r[0]=0.;
         r[nc-1]=0.;

      // update
         dt = cfl/lmax;
         for( j=0; j<nc; j++ )
        {
            qc = IdealGas2D::ConservedVariables( gas, q[j] );
            qc+= dt*r[j];
            q[j] = SolutionVariables( gas, qc );
            r[j] = 0.;
        }
     }

   // save solution
      std::ofstream solutionFile( "sods.dat" );
      IdealGas2D::State s;
      if( solutionFile.is_open() )
     {
         for( i=0; i<nc; i++ )
        {
            s = IdealGas2D::State( gas, q[i] );
            solutionFile << s.density()     << " "
                         << s.velocityX()   << " "
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

      delete[] q;
      delete[] r;

      return 0;
  }
