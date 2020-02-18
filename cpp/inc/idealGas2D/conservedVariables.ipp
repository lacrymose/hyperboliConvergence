
namespace IdealGas2D
{
   inline ConservedVariables::ConservedVariables()
  {
      var[0]=0.;
      var[1]=0.;
      var[2]=0.;
      var[3]=0.;
      return;
  }

   inline ConservedVariables ConservedVariables::operator+( ConservedVariables q )
  {
      ConservedVariables q1;
      q1[0] = var[0] + q[0];
      q1[1] = var[1] + q[1];
      q1[2] = var[2] + q[2];
      q1[3] = var[3] + q[3];
      return q1;
  }

   inline ConservedVariables ConservedVariables::operator-( ConservedVariables q )
  {
      ConservedVariables q1;
      q1[0] = var[0] - q[0];
      q1[1] = var[1] - q[1];
      q1[2] = var[2] - q[2];
      q1[3] = var[3] - q[3];
      return q1;
  }

   inline ConservedVariables& ConservedVariables::operator=( float d )
  {
      var[0]=d;
      var[1]=d;
      var[2]=d;
      var[3]=d;
      return *this;
  }

   inline ConservedVariables& ConservedVariables::operator+=( ConservedVariables q )
  {
      var[0]+=q[0];
      var[1]+=q[1];
      var[2]+=q[2];
      var[3]+=q[3];
      return *this;
  }

   inline ConservedVariables& ConservedVariables::operator-=( ConservedVariables q )
  {
      var[0]-=q[0];
      var[1]-=q[1];
      var[2]-=q[2];
      var[3]-=q[3];
      return *this;
  }

   inline ConservedVariables& ConservedVariables::operator*=( float d )
  {
      var[0]*=d;
      var[1]*=d;
      var[2]*=d;
      var[3]*=d;
      return *this;
  }

   inline ConservedVariables& ConservedVariables::operator/=( float d )
  {
      float d1=1./d;
      var[0]*=d1;
      var[1]*=d1;
      var[2]*=d1;
      var[3]*=d1;
      return *this;
  }

   inline ConservedVariables operator*( ConservedVariables q0, float d )
  {
      ConservedVariables q1;
      q1[0] = d*q0[0];
      q1[1] = d*q0[1];
      q1[2] = d*q0[2];
      q1[3] = d*q0[3];
      return q1;
  }

   inline ConservedVariables operator*( float d, ConservedVariables q0 )
  {
      ConservedVariables q1;
      q1[0] = d*q0[0];
      q1[1] = d*q0[1];
      q1[2] = d*q0[2];
      q1[3] = d*q0[3];
      return q1;
  }

   inline std::ostream& operator<<( std::ostream& os, ConservedVariables q )
  {
      os << q[0] << " "
         << q[1] << " "
         << q[2] << " "
         << q[3];
      return os;
  }
}

