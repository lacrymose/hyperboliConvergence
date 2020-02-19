# include <idealGas2D/idealGas2D.h>

# include <iostream>
# include <assert.h>
# include <math.h>

using namespace IdealGas2D;
   int main()
  {
      ConservedVariables qc0,qc1;
      ViscousVariables   qv0,qv1;

      ViscousVariables   dqv;

      Species gas;
      gas.air();

      qc0 = 1.8;
      for( int i=0; i<4; i++ )
     {
//       std::cout << qc0[i] << std::endl;
         assert( fabs( qc0[i] - 1.8 ) < 1e-6 );
     }

      qc0[0]=11.;
      qc0[1]=12.;
      qc0[2]=13.;
      qc0[3]=14.;

      qv0 = ViscousVariables(   gas, qc0 );
      qc1 = ConservedVariables( gas, qv0 );

      for( int i=0; i<4; i++ )
     {
         assert( fabs( qc0[i] - qc1[i] ) < 1e-6 );
     }

      IdealGas2D::State state_c0( gas, qc0 );
      IdealGas2D::State state_v0( gas, qv0 );

      qc1 = ConservedVariables( gas, state_c0 );
      qv1 = ViscousVariables(   gas, state_v0 );

      for( int i=0; i<4; i++ )
     {
         assert( fabs( qc0[i] - qc1[i] ) < 1e-6 );
         assert( fabs( qv0[i] - qv1[i] ) < 1e-6 );
     }

      for( int i=0; i<8; i++ )
     {
//       std::cout << state_c0.state[i] << "   " << state_v0.state[i] << std::endl;

         assert( fabs( state_c0.state[i] - state_v0.state[i] ) < 1e-6 );
     }

      return 0;
  }
