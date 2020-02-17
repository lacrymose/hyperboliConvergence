# include "cppunit/ui/text/TestRunner.h"
# include "cppunit/TestResult.h"

# include "test-BLANK.h"

   int main()
  {
      CppUnit::TextUi::TestRunner   runner;

      runner.addTest( Test_BLANK::suite() );

      bool wasSuccessful = runner.run( "", false );

      return !wasSuccessful;
  }
