# ifndef  TEST_BLANK_H
# define  TEST_BLANK_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include <BLANK.h>

/*
   Tests for public methods of the BLANK class
*/

   class Test_BLANK : public CppUnit::TestFixture
  {
   private:
      CPPUNIT_TEST_SUITE( Test_BLANK );

         CPPUNIT_TEST( test_BLANK );

      CPPUNIT_TEST_SUITE_END();

   public:
      void setUp();
      void tearDown();

      void test_BLANK();
  };

CPPUNIT_TEST_SUITE_REGISTRATION( Test_BLANK );

# endif
