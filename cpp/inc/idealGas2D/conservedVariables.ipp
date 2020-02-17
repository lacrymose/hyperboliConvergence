
namespace IdealGas2D
{
   inline ConservedVariables ConservedVariables::operator+( ConservedVariables& q )
  {
      ConservedVariables q1;
      q1[0] = var[0] + q[0];
      q1[1] = var[1] + q[1];
      q1[2] = var[2] + q[2];
      q1[3] = var[3] + q[3];
      return q1;
  }

   inline ConservedVariables ConservedVariables::operator-( ConservedVariables& q )
  {
      ConservedVariables q1;
      q1[0] = var[0] - q[0];
      q1[1] = var[1] - q[1];
      q1[2] = var[2] - q[2];
      q1[3] = var[3] - q[3];
      return q1;
  }

   inline ConservedVariables& ConservedVariables::operator+=( ConservedVariables& q )
  {
      var[0]+=q[0];
      var[1]+=q[1];
      var[2]+=q[2];
      var[3]+=q[3];
      return *this;
  }

   inline ConservedVariables& ConservedVariables::operator-=( ConservedVariables& q )
  {
      var[0]-=q[0];
      var[1]-=q[1];
      var[2]-=q[2];
      var[3]-=q[3];
      return *this;
  }

   inline std::ostream& operator<<( std::ostream& os, ConservedVariables& q )
  {
      os << q[0] << " "
         << q[1] << " "
         << q[2] << " "
         << q[3] << std::endl;
      return os;
  }
}

