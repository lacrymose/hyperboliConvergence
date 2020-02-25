
# include <idealGas2D/idealGas2D.h>
# include <array1D/array1D.h>

# include <iostream>
# include <assert.h>
# include <math.h>

   int main()
  {
      int i,j;
      IdealGas2D::Species  gas;
      IdealGas2D::VariableSet<IdealGas2D::VariableType<'c'>>     qc0,qc1,qc2;
      IdealGas2D::VariableDelta<IdealGas2D::VariableType<'c'>>  dqc0,dqc1;

      gas.air();

      qc0=1.8;
//    std::cout << qc0 << std::endl;
      for( i=0; i<4; i++ )
     {
         assert( fabs( qc0[i]-1.8 ) < 1e-6 );
     }

      dqc0=1.8;
      dqc0/=1.35;
//    std::cout << dqc0 << std::endl;
      for( i=0; i<4; i++ )
     {
         assert( fabs( dqc0[i]-1.333333333 ) < 1e-6 );
     }

      dqc0*=0.75;
//    std::cout << dqc0 << std::endl;
      for( i=0; i<4; i++ )
     {
         assert( fabs( dqc0[i]-1. ) < 1e-6 );
     }

      dqc1=3.2;
      qc0+=dqc1;
//    std::cout << qc0 << std::endl;
      for( i=0; i<4; i++ )
     {
         assert( fabs( qc0[i]-5. ) < 1e-6 );
     }

      qc0-=dqc0;
//    std::cout << qc0 << std::endl;
      for( i=0; i<4; i++ )
     {
         assert( fabs( qc0[i]-4. ) < 1e-6 );
     }

      qc1 = 0.83;
      dqc1 = dqc0*2;
//    std::cout << dqc0 << std::endl;
//    std::cout << dqc1 << std::endl;
      for( i=0; i<4; i++ )
     {
         assert( fabs( dqc1[i]-2. ) < 1e-6 );
     }

      dqc1 = 2*dqc0;
//    std::cout << dqc0 << std::endl;
//    std::cout << dqc1 << std::endl;
      for( i=0; i<4; i++ )
     {
         assert( fabs( dqc1[i]-2. ) < 1e-6 );
     }

      dqc1 = dqc0/2;
//    std::cout << dqc1 << std::endl;
      for( i=0; i<4; i++ )
     {
         assert( fabs( dqc1[i]-0.5 ) < 1e-6 );
     }

//    std::cout << std::endl;

      int n=3;
      Array::Array1D<IdealGas2D::ConservedVariables> arrayCV(n);
      Array::Array1D<IdealGas2D::ConservedVariables> arrayCW(n);
      Array::Array1D<IdealGas2D::ConservedDelta>    arrayDCV(n);

      arrayCV = 6.7;
      for( j=0; j<n; j++ )
     {
//       std::cout << arrayCV[j] << std::endl;
         for( i=0; i<4; i++ )
        {
            assert( fabs( arrayCV[j][i]-6.7 ) < 1e-6 );
        }
     }
//    std::cout << std::endl;

      arrayDCV = 3.3;
      arrayCV+=arrayDCV;
      for( j=0; j<n; j++ )
     {
//       std::cout << arrayCV[j] << std::endl;
         for( i=0; i<4; i++ )
        {
            assert( fabs( arrayCV[j][i]-10. ) < 1e-6 );
        }
     }
//    std::cout << std::endl;

      arrayCV-=arrayDCV;
      for( j=0; j<n; j++ )
     {
//       std::cout << arrayCV[j] << std::endl;
         for( i=0; i<4; i++ )
        {
            assert( fabs( arrayCV[j][i]-6.7 ) < 1e-6 );
        }
     }
//    std::cout << std::endl;

      arrayCW=2.1;
      arrayDCV=arrayCV-arrayCW;
      for( j=0; j<n; j++ )
     {
//       std::cout << arrayCV[j] << std::endl;
         for( i=0; i<4; i++ )
        {
            assert( fabs( arrayDCV[j][i]-4.6 ) < 1e-6 );
        }
     }
//    std::cout << std::endl;

      arrayCV+=arrayDCV;
      for( j=0; j<n; j++ )
     {
//       std::cout << arrayCV[j] << std::endl;
         for( i=0; i<4; i++ )
        {
            assert( fabs( arrayCV[j][i]-11.3 ) < 1e-6 );
        }
     }
//    std::cout << std::endl;

      arrayCV-=arrayDCV;
      for( j=0; j<n; j++ )
     {
//       std::cout << arrayCV[j] << std::endl;
         for( i=0; i<4; i++ )
        {
            assert( fabs( arrayCV[j][i]-6.7 ) < 1e-6 );
        }
     }
//    std::cout << std::endl;

      arrayDCV=arrayCV-arrayCW;
      for( j=0; j<n; j++ )
     {
//       std::cout << arrayCV[j] << std::endl;
         for( i=0; i<4; i++ )
        {
            assert( fabs( arrayDCV[j][i]-4.6 ) < 1e-6 );
        }
     }
//    std::cout << std::endl;

      arrayCV=arrayCV+arrayDCV;
      for( j=0; j<n; j++ )
     {
//       std::cout << arrayCV[j] << std::endl;
         for( i=0; i<4; i++ )
        {
            assert( fabs( arrayCV[j][i]-11.3 ) < 1e-6 );
        }
     }
//    std::cout << std::endl;

      arrayCV=arrayCV-arrayDCV;
      for( j=0; j<n; j++ )
     {
//       std::cout << arrayCV[j] << std::endl;
         for( i=0; i<4; i++ )
        {
            assert( fabs( arrayCV[j][i]-6.7 ) < 1e-6 );
        }
     }
//    std::cout << std::endl;

      return 0;
  }
