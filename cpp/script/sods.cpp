# include <idealGas2D/idealGas2D.h>

# include <iostream>
# include <fstream>
# include <assert.h>
# include <math.h>

   int main()
  {
      std::ifstream initialStatesFile( "initial.dat" );
      std::ifstream parametersFile( "parameters.dat" );

      IdealGas2D::ViscousVariables ql,qr;   // left/right initial states
      int      nc,nt;   // number of cells and number of timesteps
      float      cfl;   // cfl number for timestepping

      int   i,j,k;
      IdealGas2D::Species  gas;

      IdealGas2D::ConservedVariables f,qc;
      IdealGas2D::ViscousVariables     qv;
      float l,lmax,dt;
      float n[3]={1,0,1};

      gas.air();

   // read parameters
      if( initialStatesFile.is_open() )
     {
         for( i=0; i<4; i++ ){ initialStatesFile >> ql[i]; }
         for( i=0; i<4; i++ ){ initialStatesFile >> qr[i]; }
     }
      else
     {
         std::cout << "cannot open \"initial.dat\" for reading\n" << std::endl;
         return 1;
     }

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
   // ensure domain is symmetric and waves will not reach boundaries
      assert( nc%2==0 );
      assert( nt*cfl<0.5*nc );

   // initialise
      // solution and residual arrays
      IdealGas2D::ViscousVariables   *q;
      IdealGas2D::ConservedVariables *r;

      q = new IdealGas2D::ViscousVariables[  nc];
      r = new IdealGas2D::ConservedVariables[nc];

      for( i=0;    i<nc/2; i++ ){ q[i]=ql; }
      for( i=nc/2; i<nc;   i++ ){ q[i]=qr; }

   // timesteps
      lmax=-1;
      l=-1;
      for( i=0; i<nt; i++ )
     {
      // flux residual

         // left boundary
         IdealGas2D::ausm( gas, n, q[0],q[0], f,l );
         lmax = fmax( lmax, l );
         r[0]+= f;

         // interior
         for( j=0; j<nc-1; j++ )
        {
            IdealGas2D::ausm( gas, n, q[j],q[j+1], f,l );
            lmax = fmax( lmax, l );
            r[j]-= f;
            r[j]+= f;
        }

         // right boundary
         IdealGas2D::ausm( gas, n, q[nc-1],q[nc-1], f,l );
         lmax = fmax( lmax, l );
         r[nc-1]-=f;

      // update
         dt = cfl/lmax;
         for( j=0; j<nc; j++ )
        {
            qc = IdealGas2D::conservedVariables( gas, q[j] );
            for( k=0; k<4; k++ ){ qc[k]+=dt*f[k]; }
            qv = IdealGas2D::viscousVariables( gas, qc );
            q[j] = qv;
        }
     }

   // save solution
      std::ofstream solutionFile( "sods.dat" );
      if( solutionFile.is_open() )
     {
         for( i=0; i<nc; i++ )
        {
            solutionFile << q[i];
        }
     }
      else
     {
         std::cout << "cannot open \"sods.dat for writing\"\n" << std::endl;
         return 1;
     }

      return 0;
  }
