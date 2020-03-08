
# include <idealGas2D/test-State.h>

# include <types.h>

   void Test_IdealGas2D_State::setUp(){}

   void Test_IdealGas2D_State::tearDown(){}

   void Test_IdealGas2D_State::test_accessors()
  {
      s0.state[0]=1.2;
      s0.state[1]=2.2;
      s0.state[2]=3.2;
      s0.state[3]=4.2;
      s0.state[4]=5.2;
      s0.state[5]=6.2;
      s0.state[6]=7.2;
      s0.state[7]=8.2;

      CPPUNIT_ASSERT_DOUBLES_EQUAL( s0.velocityX(),             1.2, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( s0.velocityY(),             2.2, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( s0.pressure(),              3.2, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( s0.density(),               4.2, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( s0.temperature(),           5.2, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( s0.specificTotalEnthalpy(), 6.2, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( s0.velocity2(),             7.2, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( s0.speedOfSound2(),         8.2, Types::EPS );
  }

   void Test_IdealGas2D_State::test_equality()
  {
      s0.state[0]=1.2;      s1.state[0]=s0.state[0];
      s0.state[1]=2.2;      s1.state[1]=s0.state[1];
      s0.state[2]=3.2;      s1.state[2]=s0.state[2];
      s0.state[3]=4.2;      s1.state[3]=s0.state[3];
      s0.state[4]=5.2;      s1.state[4]=s0.state[4];
      s0.state[5]=6.2;      s1.state[5]=s0.state[5];
      s0.state[6]=7.2;      s1.state[6]=s0.state[6];
      s0.state[7]=8.2;      s1.state[7]=s0.state[7];

   // equal states
      CPPUNIT_ASSERT(   s0==s1  );
      CPPUNIT_ASSERT( !(s0!=s1) );

   // unequal states
      s0.state[0]=0.;
      CPPUNIT_ASSERT( !(s0==s1) );
      CPPUNIT_ASSERT(   s0!=s1  );
      s0.state[0]=s1.state[0];

      s1.state[1]=1.;
      CPPUNIT_ASSERT( !(s0==s1) );
      CPPUNIT_ASSERT(   s0!=s1  );
      s0.state[1]=s1.state[1];

      s0.state[2]=2.;
      CPPUNIT_ASSERT( !(s0==s1) );
      CPPUNIT_ASSERT(   s0!=s1  );
      s0.state[2]=s1.state[2];

      s0.state[3]=3.;
      CPPUNIT_ASSERT( !(s0==s1) );
      CPPUNIT_ASSERT(   s0!=s1  );
      s0.state[3]=s1.state[3];

      s0.state[4]=4.;
      CPPUNIT_ASSERT( !(s0==s1) );
      CPPUNIT_ASSERT(   s0!=s1  );
      s0.state[4]=s1.state[4];

      s0.state[5]=5.;
      CPPUNIT_ASSERT( !(s0==s1) );
      CPPUNIT_ASSERT(   s0!=s1  );
      s0.state[5]=s1.state[5];

      s0.state[6]=6.;
      CPPUNIT_ASSERT( !(s0==s1) );
      CPPUNIT_ASSERT(   s0!=s1  );
      s0.state[6]=s1.state[6];

      s0.state[7]=7.;
      CPPUNIT_ASSERT( !(s0==s1) );
      CPPUNIT_ASSERT(   s0!=s1  );
      s0.state[7]=s1.state[7];

  }

   void Test_IdealGas2D_State::test_initialisation_conserved()
  {
  }

   void Test_IdealGas2D_State::test_initialisation_viscous()
  {
  }

