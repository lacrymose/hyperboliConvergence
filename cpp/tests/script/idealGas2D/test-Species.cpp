# include "cppunit/ui/text/TestRunner.h"
# include "cppunit/TestResult.h"

# include "idealGas2D/test-Species.h"

   int main()
  {
      CppUnit::TextUi::TestRunner   runner;

      runner.addTest( Test_IdealGas2D_Species::suite() );

      bool wasSuccessful = runner.run( "", false );

      return !wasSuccessful;
  }
