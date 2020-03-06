
# include <fieldOperations/fieldOperations.h>

# include <array1D/array1D.h>
# include <controls.h>

# include <types.h>

# include <iostream>
# include <fstream>
# include <cmath>

   int main()
  {
   // set up parameters
      int  n;  // number of points
      int  omega;  // wavenumber

      Types::Real rn0,rn1;

      Types::Real   alpha,beta,gamma;
      Types::Real  ampe,ampi,ampei,h;

      Array::Array1D<Types::Real>   r0;
      Array::Array1D<Types::Real>   re;
      Array::Array1D<Types::Real>   ri;
      Array::Array1D<Types::Real>   rei;

      Controls::GridControls1D    grid;

   // initialise grid and solution
      std::ifstream  inputFile( "data/smooth/parameters.dat" );
      if( inputFile.is_open() )
     {
         inputFile >> n;
         inputFile >> omega;
         inputFile >> alpha;
     }
      else
     {
         std::cout << "cannot open \"parameters.dat\" for reading\n" << std::endl;
         return 1;
     }
      inputFile.close();

   // initialise fourier mode
      grid.n=n;
      grid.boundaryCondition='p';
      h = 2.*M_PI/n;

      r0.resize(n);
      re.resize(n);
      ri.resize(n);
      rei.resize(n);

   // smoothing parameters
      beta = 0.25*( alpha*alpha - 1. );
      gamma = 0.25*( (1.+4.*beta)/alpha - 1. );

   // results file
      std::ofstream outputFile( "data/smooth/smooth.dat" );
      if( !outputFile.is_open() )
     {
         std::cout << "cannot open \"smooth.dat\" for writing\n" << std::endl;
         return 1;
     }

      rn0=1.0;
      int j=1;
      for( j=1; j<n/2; j++ )
     {
         omega = j;

         ampe= 1. + 4.*gamma*sin(0.5*omega*h)*sin(0.5*omega*h);
         ampi= 1. + 4.*beta*sin( 0.5*omega*h)*sin(0.5*omega*h);
         ampi= alpha/ampi;
         ampei=ampe*ampi;
         outputFile << ampe << " " << ampi << " " << ampei;

      // smooth solution
         Types::Real x;

         for( int i=0; i<n; i++ )
        {
            x = i/Types::Real(n);
            r0[i]=rn0*std::sin( 2*M_PI*omega*x );
        }

         centralExplicitSmoothing( grid,alpha, r0, re  );
         centralImplicitSmoothing( grid,alpha, r0, ri  );
         centralImplicitSmoothing( grid,alpha, re, rei );

         rn1=0.;
         for( int i=0; i<n; i++ ){ rn1=fmax(rn1,fabs(re[i])); }
         outputFile << " " << rn1/rn0;
         rn1=0.;
         for( int i=0; i<n; i++ ){ rn1=fmax(rn1,fabs(ri[i])); }
         outputFile << " " << rn1/rn0;
         rn1=0.;
         for( int i=0; i<n; i++ ){ rn1=fmax(rn1,fabs(rei[i])); }
         outputFile << " " << rn1/rn0;
         outputFile << std::endl;

     }

   // print output
//    for( int i=0; i<n; i++ )
//   {
//       outputFile << r0[i] << " " << ri[i] << std::endl;
//   }

      outputFile.close();
      return 0;
  }
