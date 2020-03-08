# ifndef  TEST_IDEALGAS2D_SPECIES_H
# define  TEST_IDEALGAS2D_SPECIES_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include <idealGas2D/idealGas2D.h>

/*
   Tests for public methods of the IdealGas2D::Species class
*/

   class Test_IdealGas2D_Species : public CppUnit::TestFixture
  {
   private:
      CPPUNIT_TEST_SUITE( Test_IdealGas2D_Species );

         CPPUNIT_TEST( test_air );

      CPPUNIT_TEST_SUITE_END();

      IdealGas2D::Species  gas;

   public:
      void setUp();
      void tearDown();

   // test initialising the species to air at standard conditions
      void test_air();
  };

CPPUNIT_TEST_SUITE_REGISTRATION( Test_IdealGas2D_Species );

# endif
