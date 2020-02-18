# include <idealGas2D/idealGas2D.h>

namespace IdealGas2D
{
   void Species::air()
  {
      pr=       0.7;
      minf=     1.0;
      gamma=    1.4;
      Rgas= 287.058;
      nu=1.81*10e-5;
  }
}
