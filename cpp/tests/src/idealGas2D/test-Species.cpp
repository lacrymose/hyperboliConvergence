
# include <idealGas2D/test-Species.h>

# include <types.h>

   void Test_IdealGas2D_Species::setUp(){}

   void Test_IdealGas2D_Species::tearDown(){}

   void Test_IdealGas2D_Species::test_air()
  {
      Types::Real prandtl=0.7;
      Types::Real machRef=1.0;
      Types::Real gamma  =1.4;
      Types::Real rGas   =287.058;
      Types::Real viscous=1.81e-5;

      gas.air();

      CPPUNIT_ASSERT_DOUBLES_EQUAL( gas.pr,    prandtl, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( gas.minf,  machRef, Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( gas.gamma, gamma,   Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( gas.Rgas,  rGas,    Types::EPS );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( gas.nu,    viscous, Types::EPS );
  }

